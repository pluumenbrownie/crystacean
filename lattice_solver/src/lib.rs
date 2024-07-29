#![warn(clippy::pedantic, clippy::nursery)]
#![allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::cast_possible_truncation,
    clippy::wildcard_imports,
    clippy::must_use_candidate,
    clippy::return_self_not_must_use
)]
// Optimisation Ideas:
// - Turn current/next generation into hashsets to remove archive requirement
//     - wont work with depth first
// - Try depth first and speed that up
// - parallellise
//     - try rayon
// - get profiler to work

use fixedbitset::FixedBitSet;
use itertools::{izip, Itertools};
use json::{object, JsonValue};
use std::{ffi::OsString, fs::File, iter::zip, sync::Arc};

use kiddo::{float::kdtree::KdTree, SquaredEuclidean};
use std::io::prelude::*;

mod points;
use points::*;

mod bit_representation;
pub use bit_representation::*;

mod inserters;
use inserters::*;

mod crown;
use crown::*;

pub mod close_vector_tree;
// use close_vector_tree::*;

pub mod site_filter;
pub mod test_points;

const BINSIZE: usize = 129;

pub struct Lattice {
    points: Vec<Arc<LatticePoint>>,
    oxygens: Vec<Oxygen>,
    source_file: Option<JsonValue>,
}

impl Lattice {
    /// Create an empty `Lattice`
    fn new() -> Self {
        Self {
            points: vec![],
            oxygens: vec![],
            source_file: None,
        }
    }

    /// Push a point to `self.points`
    fn add_point(&mut self, new_point: Arc<LatticePoint>) {
        self.points.push(new_point);
    }

    fn generate_exclusions(&mut self) {
        let oxygen_amount = self.oxygens.len();
        let mut distances_matrix = vec![];

        for number in 0..oxygen_amount {
            let mut distances_row = vec![];
            for other in 0..oxygen_amount {
                distances_row.push(self.distance_between(OxygenIndex(number), OxygenIndex(other)));
            }
            distances_matrix.push(distances_row);
        }

        for (number, oxygen) in self.oxygens.iter().enumerate() {
            for index in oxygen.sitetype.iter() {
                let mut connections = self.points[index.0].get_connections().write().unwrap();
                connections.push(OxygenIndex(number));
            }
        }
        for (number, oxygen) in self.oxygens.iter_mut().enumerate() {
            for connection in oxygen.sitetype.iter() {
                let point = &self.points[connection.0];

                #[allow(clippy::significant_drop_in_scrutinee)]
                for ox_con in point.get_connections().read().unwrap().iter() {
                    oxygen.exclusions.push(*ox_con);
                }
            }
            // points which are to close to one another should exclude eachother
            for other in 0..oxygen_amount {
                if distances_matrix[number][other] < 0.02 && number != other {
                    oxygen.exclusions.push(OxygenIndex(other));
                };
            }
            oxygen.exclusions.dedup();
        }
    }

    /// Turns the lattice problem in an abstracted form based on bitarrays.
    fn generate_intermediary(&self, options: BitArraySettings) -> BitArrayRepresentation {
        let filled_sites = FixedBitSet::with_capacity(self.oxygens.len());
        let mut exclusion_matrix: Vec<FixedBitSet> = vec![];
        let mut distances_matrix = vec![];

        let mut tripoint_mask = FixedBitSet::with_capacity(self.oxygens.len());
        let mut midpoint_mask = FixedBitSet::with_capacity(self.oxygens.len());
        let mut singlet_mask = FixedBitSet::with_capacity(self.oxygens.len());

        for (number, oxygen) in self.oxygens.iter().enumerate() {
            let mut exclusions = FixedBitSet::with_capacity(self.oxygens.len());
            for exclusion in &oxygen.exclusions {
                exclusions.set(exclusion.0, true);
            }

            let mut distances_row = vec![];
            for other in 0..self.oxygens.len() {
                distances_row.push(if other > number {
                    self.distance_between(OxygenIndex(number), OxygenIndex(other))
                } else {
                    0.0
                });
            }
            distances_matrix.push(distances_row);

            match &oxygen.sitetype {
                SiteType::Tripoint(_) => tripoint_mask.set(number, true),
                SiteType::Midpoint(_) => midpoint_mask.set(number, true),
                SiteType::Singlet(_) => singlet_mask.set(number, true),
            };
            exclusion_matrix.push(exclusions);
        }

        BitArrayRepresentation {
            filled_sites,
            exclusion_matrix,
            distances_matrix,
            tripoint_mask,
            midpoint_mask,
            singlet_mask,
            filter: None,
            options,
        }
    }

    pub fn find_max(&self) -> (f32, f32) {
        let (mut max_x, mut max_y) = (0.0, 0.0);

        for point in &self.points {
            if point.x > max_x && point.y == 0.0 {
                max_x = point.x;
            }
            if point.y > max_y {
                max_y = point.y;
            }
        }

        (max_x, max_y)
    }

    pub fn distance_between(&self, index_one: OxygenIndex, index_two: OxygenIndex) -> f32 {
        let (max_x, max_y) = self.find_max();
        let one = &self.oxygens[index_one.0];
        let two = &self.oxygens[index_two.0];

        let delta_x = {
            let dx = (one.x - two.x).abs();
            if dx > max_x / 2.0 {
                dx - max_x
            } else {
                dx
            }
        };
        let delta_y = {
            let dy = (one.y - two.y).abs();
            if dy > max_y / 2.0 {
                dy - max_y
            } else {
                dy
            }
        };

        #[allow(clippy::suboptimal_flops)]
        (delta_x.powi(2) + delta_y.powi(2) + (one.z - two.z).powi(2)).sqrt()
    }

    /// `distance_margin` should be 1.1 for 2D, 1.4 for 3D
    #[must_use]
    pub fn python_new(
        input_lattice: Vec<(Vec<f32>, Vec<Vec<f32>>)>,
        distance_margin: f32,
        autodetect_margin: bool,
    ) -> Self {
        // Convert 2D structures to 3D
        let lattice_3d = turn_2d_3d(input_lattice);

        // Create the silicon lattice
        let mut out_lattice = create_silicon_lattice(lattice_3d);

        let first_point_location = {
            let first_point = &out_lattice.points[0];
            [first_point.x, first_point.y, first_point.z]
        };

        // Fill in the oxygens
        let silicon_iterator = out_lattice
            .points
            .iter()
            .map(|p| [p.x, p.y, p.z])
            .collect_vec();

        let kdtree: KdTree<_, u64, 3, BINSIZE, u32> = KdTree::from(&silicon_iterator);
        let node_search_distance = if autodetect_margin {
            kdtree.nearest_n::<SquaredEuclidean>(&first_point_location, 2)[1].distance
                * distance_margin
        } else {
            distance_margin.powi(2)
        };

        insert_tripoints(&mut out_lattice, &kdtree, node_search_distance);
        insert_midpoints(&mut out_lattice, &kdtree, node_search_distance);
        insert_singles(&mut out_lattice);

        out_lattice.generate_exclusions();
        out_lattice
    }

    fn add_source_file(&mut self, source_file: JsonValue) {
        self.source_file = Some(source_file);
    }

    /// Create a `Lattice` from an ASE json file.
    ///
    /// Usefull for creating secondary layers on top of structures processed with CP2K.
    ///
    /// # Panics
    /// Will panic when supplied with missing or invalid file.
    #[must_use]
    pub fn from_dft_json(filename: String, distance_margin: f32, autodetect_margin: bool) -> Self {
        let mut buffer = String::new();
        let mut file = File::open(filename).expect("Opening file failed.");
        file.read_to_string(&mut buffer)
            .expect("Reading file failed.");
        let parsed = json::parse(&buffer).expect("Parsing file failed.");

        let last_id = &parsed["ids"][parsed["ids"].len() - 1].to_string();

        // let hydrogen_amount = &parsed[last_id]["numbers"]["__ndarray__"][2]
        //     .members()
        //     .filter_map(|n| n.as_usize())
        //     .filter(|n| n == &1)
        //     .count();

        let numbers = &parsed[last_id]["positions"]["__ndarray__"][2];
        let atoms = numbers
            .members()
            .map(|j| j.as_f32().unwrap())
            .collect_vec()
            .chunks_exact(3)
            .map(<[f32]>::to_vec)
            .zip(
                parsed[last_id]["numbers"]["__ndarray__"][2]
                    .members()
                    .map(json::JsonValue::as_i32),
            )
            .collect_vec();

        let hydrogenated_ends = atoms
            .iter()
            .filter(|(_, i)| i == &Some(1))
            .map(|(v, _)| v)
            .filter(|v| v[2] < 20.0)
            .collect_vec();

        let cell = &parsed[last_id]["cell"]["array"]["__ndarray__"][2];
        let (x_vec, y_vec, _z_vec) = cell
            .members()
            .map(|v| v.as_f32().unwrap())
            .tuples::<(_, _, _)>()
            .collect_tuple()
            .unwrap_or_else(|| panic!("Json 'cell' property format is incorrect: {cell:?}"));

        let input_lattice = {
            let mut input_lattice = vec![];
            for end in hydrogenated_ends {
                let mut new_point = (end.clone(), vec![]);
                new_point
                    .1
                    .push(vec![end[0] + x_vec.0, end[1] + x_vec.1, end[2] + x_vec.2]);
                new_point
                    .1
                    .push(vec![end[0] + y_vec.0, end[1] + y_vec.1, end[2] + y_vec.2]);
                new_point.1.push(vec![
                    end[0] + x_vec.0 + y_vec.0,
                    end[1] + x_vec.1 + y_vec.1,
                    end[2] + x_vec.2 + y_vec.2,
                ]);
                new_point.1.push(vec![
                    end[0] + x_vec.0 - y_vec.0,
                    end[1] + x_vec.1 - y_vec.1,
                    end[2] + x_vec.2 - y_vec.2,
                ]);
                new_point.1.push(vec![
                    end[0] - x_vec.0 + y_vec.0,
                    end[1] - x_vec.1 + y_vec.1,
                    end[2] - x_vec.2 + y_vec.2,
                ]);
                input_lattice.push(new_point);
            }
            input_lattice
        };

        let mut lattice = Self::python_new(input_lattice, distance_margin, autodetect_margin);
        lattice.add_source_file(parsed.clone());

        lattice
    }

    /// Output the `Lattice` for diagnostic purpouses.
    ///
    /// # Panics
    /// Can panic when `Lattice` was not created from file or file writing fails.
    pub fn diagnostic_ase(&self) {
        let parsed = self.source_file.as_ref().unwrap();
        let oxygens = &self.oxygens;
        let last_id = &parsed["ids"][parsed["ids"].len() - 1].to_string();

        let mut new_numbers = parsed[last_id]["numbers"].clone();
        let mut new_positions = parsed[last_id]["positions"].clone();

        println!("Added {} oxygens.", oxygens.len());

        for oxygen in oxygens {
            new_numbers["__ndarray__"][2]
                .push(8)
                .expect("new_numbers[\"__ndarray__\"][2].push(8)");
            new_positions["__ndarray__"][2]
                .push(oxygen.x)
                .expect("new_positions[\"__ndarray__\"][2].push(oxygen.x)");
            new_positions["__ndarray__"][2]
                .push(oxygen.y)
                .expect("new_positions[\"__ndarray__\"][2].push(oxygen.y)");
            new_positions["__ndarray__"][2]
                .push(oxygen.z)
                .expect("new_positions[\"__ndarray__\"][2].push(oxygen.z)");
        }
        new_numbers["__ndarray__"][0][0] = new_numbers["__ndarray__"][2].len().into();
        new_positions["__ndarray__"][0][0] = new_numbers["__ndarray__"][2].len().into();

        let mut export_data = json::JsonValue::new_object();
        export_data["1"] = object! {
            cell: parsed[last_id]["cell"].clone(),
            ctime: parsed[last_id]["ctime"].clone(),
            mtime: parsed[last_id]["mtime"].clone(),
            numbers: new_numbers,
            pbc: parsed[last_id]["pbc"].clone(),
            positions: new_positions,
            unique_id: "Not unique",
            user: parsed[last_id]["user"].clone(),
        };

        let mut file = File::create("output.json").unwrap();
        file.write_all(export_data.pretty(4).as_bytes()).unwrap();
    }

    /// Helper method for the `no_rings` methods.
    fn no_rings_detector(&self) -> (Vec<[f32; 3]>, Vec<usize>) {
        let parsed = self.source_file.as_ref().unwrap();

        let last_id = &parsed["ids"][parsed["ids"].len() - 1].to_string();

        let numbers = &parsed[last_id]["positions"]["__ndarray__"][2];
        let elements = &parsed[last_id]["numbers"]["__ndarray__"][2];
        let atoms = numbers
            .members()
            .map(|j| j.as_f32().unwrap())
            .collect_vec()
            .chunks_exact(3)
            .map(<[f32]>::to_vec)
            .collect_vec();

        let top_silicon_locations = atoms
            .iter()
            .zip_eq(elements.members())
            .sorted_by(|&a, &b| a.0[2].total_cmp(&b.0[2]))
            .filter(|(s, e)| e.as_i64().unwrap() == 14 && (s[2] < 10.0))
            .map(|(v, _)| [v[0], v[1], v[2]])
            .collect_vec();

        let points_vector = self.points.iter().map(|p| [p.x, p.y, p.z]).collect_vec();
        let points_tree: KdTree<_, u64, 3, BINSIZE, u32> = (&points_vector).into();

        let point_group_vector = {
            let mut point_group_vector = vec![1000; self.points.len()];

            for (number, atom) in top_silicon_locations.iter().enumerate() {
                let close_points = points_tree.within::<SquaredEuclidean>(atom, 1.7f32.powi(2));
                // println!("{number}:{atom:?} --- {close_points:?}");
                for point in close_points {
                    point_group_vector[point.item as usize] = number;
                }
            }
            point_group_vector
        };
        (points_vector, point_group_vector)
    }

    /// Returns a numbered list of ring locations.
    ///
    /// # Panics
    /// Can panic when `Lattice` was not created from file.
    #[must_use]
    pub fn no_rings_plot(&self) -> Vec<(usize, f32, f32)> {
        let (points_vector, point_group_vector) = self.no_rings_detector();
        println!("points_vector = {points_vector:?}");
        println!("length: {}", points_vector.len());
        println!("point_group_vector = {point_group_vector:?}");

        izip!(
            point_group_vector,
            self.points.iter().map(|o| o.x),
            self.points.iter().map(|o| o.y)
        )
        .collect_vec()
    }

    /// Create a `SiteFilter` to prevent small rings from forming in the material.
    ///
    /// # Panics
    /// Can panic when `Lattice` was not created from file.
    #[must_use]
    pub fn no_rings(&self) -> site_filter::SiteFilter {
        let (_, point_group_vector) = self.no_rings_detector();

        let mut disabled_oxygens = vec![];

        for (number, oxygen) in self.oxygens.iter().enumerate() {
            match oxygen.sitetype {
                SiteType::Singlet(_) => {}
                SiteType::Midpoint(p) => {
                    let connections = p.0;
                    let same_group = |a: usize, b: usize| {
                        point_group_vector[connections[a].0] == point_group_vector[connections[b].0]
                    };

                    if same_group(0, 1) {
                        disabled_oxygens.push(OxygenIndex(number));
                    };
                }
                SiteType::Tripoint(p) => {
                    let connections = p.0;
                    let same_group = |a: usize, b: usize| {
                        point_group_vector[connections[a].0] == point_group_vector[connections[b].0]
                    };

                    if same_group(0, 1) | same_group(1, 2) | same_group(0, 2) {
                        disabled_oxygens.push(OxygenIndex(number));
                    };
                }
            }
        }

        site_filter::SiteFilter {
            wrapped: disabled_oxygens,
        }
    }

    /// Returns the coordinates of the lattice points in two lists.
    /// Use with the * star operator in a `plt.plot` function:
    ///
    /// ```python
    /// plt.plot(*solved_lattice.points_to_plot(), "o")
    /// ```
    #[must_use]
    pub fn points_to_plot(&self) -> (Vec<f32>, Vec<f32>) {
        let x_points = self.points.iter().map(|p| p.x).collect_vec();
        let y_points = self.points.iter().map(|p| p.y).collect_vec();
        (x_points, y_points)
    }

    /// Returns the coordinates of the oxygen points in two lists.
    /// Use with the * star operator in a `plt.plot` function:
    ///
    /// ```python
    /// plt.plot(*solved_lattice.oxygens_to_plot(), "o")
    /// ```
    #[must_use]
    pub fn oxygens_to_plot(&self) -> (Vec<f32>, Vec<f32>) {
        let x_points = self.oxygens.iter().map(|p| p.x).collect_vec();
        let y_points = self.oxygens.iter().map(|p| p.y).collect_vec();
        (x_points, y_points)
    }

    /// Returns the coordinates of the tripoints in two lists.
    /// Use with the * star operator in a `plt.plot` function:
    ///
    /// ```python
    /// plt.plot(*solved_lattice.tripoints_to_plot(), "o")
    /// ```
    #[must_use]
    pub fn tripoints_to_plot(&self) -> (Vec<f32>, Vec<f32>) {
        let x_points = self
            .oxygens
            .iter()
            .filter(|p| matches!(p.sitetype, SiteType::Tripoint(_)))
            .map(|p| p.x)
            .collect_vec();
        let y_points = self
            .oxygens
            .iter()
            .filter(|p| matches!(p.sitetype, SiteType::Tripoint(_)))
            .map(|p| p.y)
            .collect_vec();
        (x_points, y_points)
    }

    /// Returns the coordinates of the tripoints in two lists.
    /// Use with the * star operator in a `plt.plot` function:
    ///
    /// ```python
    /// plt.plot(*solved_lattice.midpoints_to_plot(), "o")
    /// ```
    #[must_use]
    pub fn midpoints_to_plot(&self) -> (Vec<f32>, Vec<f32>) {
        let x_points = self
            .oxygens
            .iter()
            .filter(|p| matches!(p.sitetype, SiteType::Midpoint(_)))
            .map(|p| p.x)
            .collect_vec();
        let y_points = self
            .oxygens
            .iter()
            .filter(|p| matches!(p.sitetype, SiteType::Midpoint(_)))
            .map(|p| p.y)
            .collect_vec();
        (x_points, y_points)
    }

    /// Returns the coordinates of the tripoints in two lists.
    /// Use with the * star operator in a `plt.plot` function:
    ///
    /// ```python
    /// plt.plot(*solved_lattice.singlets_to_plot(), "o")
    /// ```
    #[must_use]
    pub fn singlets_to_plot(&self) -> (Vec<f32>, Vec<f32>) {
        let x_points = self
            .oxygens
            .iter()
            .filter(|p| matches!(p.sitetype, SiteType::Singlet(_)))
            .map(|p| p.x)
            .collect_vec();
        let y_points = self
            .oxygens
            .iter()
            .filter(|p| matches!(p.sitetype, SiteType::Singlet(_)))
            .map(|p| p.y)
            .collect_vec();
        (x_points, y_points)
    }

    /// Generates a more efficient representation of the lattice
    /// problem for the given lattice.
    #[must_use]
    pub fn get_intermediary(&self, options: BitArraySettings) -> BitArrayRepresentation {
        self.generate_intermediary(options)
    }

    /// Returns a solved version of the lattice. Usefull for plotting
    /// and exporting.
    #[must_use]
    pub fn to_solved_lattice(&self, solution: &BitArraySolution) -> Self {
        let mut solved_oxygens = vec![];
        for number in 0..self.oxygens.len() {
            if solution.0[number] {
                solved_oxygens.push(self.oxygens[number].clone());
            }
        }
        Self {
            points: self.points.clone(),
            oxygens: solved_oxygens,
            source_file: self.source_file.clone(),
        }
    }

    /// Export the solved lattice to a json file.
    ///
    /// - path must be a valid path name.
    /// - name should end with ".json".
    ///
    /// # Panics
    /// Can panic when path is invalid.
    pub fn export(&self, path: &OsString, name: String) {
        let mut data = json::JsonValue::new_object();
        {
            let mut points = vec![];
            for point in &self.points {
                let mut new_obj = json::JsonValue::new_object();
                new_obj["x"] = point.x.into();
                new_obj["y"] = point.y.into();
                new_obj["ghost"] = point.ghost_to.is_some().into();
                points.push(new_obj);
            }
            data["lattice_points"] = points.into();
        }

        {
            let mut tripoints = vec![];
            let mut midpoints = vec![];
            let mut singles = vec![];
            for oxygen in &self.oxygens {
                let mut new_obj = json::JsonValue::new_object();
                new_obj["x"] = oxygen.x.into();
                new_obj["y"] = oxygen.y.into();

                match oxygen.sitetype {
                    SiteType::Tripoint(_) => tripoints.push(new_obj),
                    SiteType::Midpoint(_) => midpoints.push(new_obj),
                    SiteType::Singlet(_) => singles.push(new_obj),
                };
            }
            data["tripoints"] = tripoints.into();
            data["midpoints"] = midpoints.into();
            data["singles"] = singles.into();
        }

        let mut filename = path.clone();
        filename.push("/");
        filename.push(name);

        let mut file = File::create(&filename).unwrap_or_else(|_| panic!("{:?}", &filename));

        file.write_all(data.pretty(4).as_bytes()).unwrap();
    }

    /// Export the found solution as a json file, which can be converted to other
    /// formats with ASE.
    ///
    /// # Panics
    /// Can panic when `Lattice` was not created from file or filename points to
    /// invalid path.
    pub fn export_as_ase_json(&self, filename: &String) {
        let parsed = self.source_file.as_ref().unwrap();
        let oxygens = &self.oxygens;
        let last_id = &parsed["ids"][parsed["ids"].len() - 1].to_string();

        let old_numbers = parsed[last_id]["numbers"].clone();
        let mut new_numbers = parsed[last_id]["numbers"].clone();
        let mut new_positions = parsed[last_id]["positions"].clone();

        let z_coords = new_positions["__ndarray__"][2].members().skip(2).step_by(3);

        new_numbers["__ndarray__"][2].clear();
        for (old_number, z) in zip(old_numbers["__ndarray__"][2].members(), z_coords) {
            if (old_number.as_usize() == Some(1usize)) & (z.as_f64() < Some(20.0)) {
                new_numbers["__ndarray__"][2]
                    .push(8)
                    .expect("Pushing new number failed.");
            } else {
                new_numbers["__ndarray__"][2]
                    .push(old_number.clone())
                    .expect("Pushihng old number failed");
            }
        }

        let mut c = (0usize, 0usize, 0usize);
        for oxygen in oxygens {
            new_numbers["__ndarray__"][2]
                .push(14)
                .expect("new_numbers[\"__ndarray__\"][2].push(8)");
            new_positions["__ndarray__"][2]
                .push(oxygen.x)
                .expect("new_positions[\"__ndarray__\"][2].push(oxygen.x)");
            new_positions["__ndarray__"][2]
                .push(oxygen.y)
                .expect("new_positions[\"__ndarray__\"][2].push(oxygen.y)");
            new_positions["__ndarray__"][2]
                .push(oxygen.z)
                .expect("new_positions[\"__ndarray__\"][2].push(oxygen.z)");
            match oxygen.sitetype {
                SiteType::Tripoint(_) => c = (c.0 + 1, c.1, c.2),
                SiteType::Midpoint(_) => c = (c.0, c.1 + 1, c.2),
                SiteType::Singlet(_)  => c = (c.0, c.1, c.2 + 1),
            };
            self.add_crown(oxygen, &mut new_numbers, &mut new_positions);
        }
        new_numbers["__ndarray__"][0][0] = new_numbers["__ndarray__"][2].len().into();
        new_positions["__ndarray__"][0][0] = new_numbers["__ndarray__"][2].len().into();

        let mut export_data = json::JsonValue::new_object();
        export_data["1"] = object! {
            cell: parsed[last_id]["cell"].clone(),
            ctime: parsed[last_id]["ctime"].clone(),
            mtime: parsed[last_id]["mtime"].clone(),
            numbers: new_numbers,
            pbc: parsed[last_id]["pbc"].clone(),
            positions: new_positions,
            unique_id: "Not unique",
            user: parsed[last_id]["user"].clone(),
            tripoints: c.0,
            midpoints: c.1,
            singlets: c.2,
        };

        let mut file = File::create(filename).expect("Folder does not exist!");
        file.write_all(export_data.pretty(4).as_bytes()).unwrap();
    }

    fn add_crown(
        &self,
        oxygen: &Oxygen,
        new_numbers: &mut JsonValue,
        new_positions: &mut JsonValue,
    ) {
        let double_crown;
        let single_crown;
        let new_crown = match oxygen.sitetype {
            SiteType::Singlet(_) => {
                single_crown = triple_crown_rotated(0.0);
                single_crown.iter()
            }
            SiteType::Midpoint(p) => {
                double_crown = double_crown_rotated(double_angle(
                    &self.points[p.0[0].0],
                    &self.points[p.0[1].0],
                ));
                double_crown.iter()
            }
            SiteType::Tripoint(_) => SINGLE_CROWN.iter(),
        };

        for (x, y, z) in new_crown.tuples() {
            new_numbers["__ndarray__"][2]
                .push(1)
                .expect("Adding to borrowed new_numbers failed");
            new_positions["__ndarray__"][2]
                .push(oxygen.x + x)
                .expect("Adding to borrowed new_positions failed");
            new_positions["__ndarray__"][2]
                .push(oxygen.y + y)
                .expect("Adding to borrowed new_positions failed");
            new_positions["__ndarray__"][2]
                .push(oxygen.z + z)
                .expect("Adding to borrowed new_positions failed");
        }
    }
}

fn create_silicon_lattice(lattice_3d: Vec<(Vec<f32>, Vec<Vec<f32>>)>) -> Lattice {
    let mut out_lattice = Lattice::new();

    for (location, ghosts) in lattice_3d {
        let new_point = LatticePoint::new(location[0], location[1], location[2], None);
        out_lattice.add_point(new_point.clone());

        for ghost in ghosts {
            out_lattice.add_point(LatticePoint::new(
                ghost[0],
                ghost[1],
                ghost[2],
                Some(new_point.clone()),
            ));
        }
    }
    out_lattice
        .points
        .sort_by_key(|p| 100.0f32.mul_add(p.x, p.y).round() as i32);
    out_lattice
}

fn turn_2d_3d(input_lattice: Vec<(Vec<f32>, Vec<Vec<f32>>)>) -> Vec<(Vec<f32>, Vec<Vec<f32>>)> {
    if input_lattice[0].0.len() == 2 {
        let mut new_lattice = vec![];
        for (point, ghosts) in input_lattice {
            let mut new_point = point.clone();
            new_point.push(0.0);
            let mut new_ghosts = ghosts.clone();
            for g in &mut new_ghosts {
                g.push(0.0);
            }
            new_lattice.push((new_point, new_ghosts));
        }
        new_lattice
    } else if input_lattice[0].0.len() == 3 {
        input_lattice
    } else {
        panic!("Input lattice layout is incorrect: points must be two or threedimentional.")
    }
}
