// Optimisation Ideas:
// - Turn current/next generation into hashsets to remove archive requirement
//     - wont work with depth first
// - Try depth first and speed that up
// - parallellise
//     - try rayon
// - get profiler to work

use fixedbitset::FixedBitSet;
use itertools::{zip_eq, Itertools};
use kdam::{tqdm, Colour, Spinner};
use std::{
    collections::HashSet,
    ffi::OsString,
    fs::File,
    io::{stderr, IsTerminal},
    mem,
    sync::{Arc, RwLock},
};

use kiddo::{KdTree, SquaredEuclidean};
use std::io::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct OxygenIndex(usize);

#[derive(Clone, Copy)]
struct LatticeIndex(usize);

struct LatticePoint {
    x: f32,
    y: f32,
    connected_to: RwLock<Vec<OxygenIndex>>,
    ghost_to: Option<Arc<LatticePoint>>,
}

impl LatticePoint {
    fn new(x: f32, y: f32, ghost_to: Option<Arc<LatticePoint>>) -> Arc<Self> {
        Arc::new(LatticePoint {
            x,
            y,
            connected_to: RwLock::new(vec![]),
            ghost_to,
        })
    }

    fn get_connections(&self) -> &RwLock<Vec<OxygenIndex>> {
        match &self.ghost_to {
            None => &self.connected_to,
            Some(point) => &point.connected_to,
        }
    }
}

#[derive(Clone)]
struct Oxygen {
    x: f32,
    y: f32,
    sitetype: SiteType,
    exclusions: Vec<OxygenIndex>,
}

impl Oxygen {
    fn new(x: f32, y: f32, sitetype: SiteType) -> Self {
        Oxygen {
            x,
            y,
            sitetype,
            exclusions: vec![],
        }
    }
}

#[derive(Clone, Copy)]
enum SiteType {
    Singlet(Singlet),
    Midpoint(Midpoint),
    Tripoint(Tripoint),
}

impl SiteType {
    fn iter(&self) -> std::slice::Iter<'_, LatticeIndex> {
        match self {
            Self::Tripoint(c) => c.0.iter(),
            Self::Midpoint(c) => c.0.iter(),
            Self::Singlet(c) => c.0.iter(),
        }
    }
}

#[derive(Clone, Copy)]
struct Singlet([LatticeIndex; 1]);

#[derive(Clone, Copy)]
struct Midpoint([LatticeIndex; 2]);

#[derive(Clone, Copy)]
struct Tripoint([LatticeIndex; 3]);

pub struct BitArraySolution(FixedBitSet);

pub struct BitArrayRepresentation {
    filled_sites: FixedBitSet,
    exclusion_matrix: Vec<FixedBitSet>,
    tripoint_mask: FixedBitSet,
    midpoint_mask: FixedBitSet,
    singlet_mask: FixedBitSet,
}

impl BitArrayRepresentation {
    fn matrix_vector_multiply(&self, vector: &FixedBitSet) -> FixedBitSet {
        let mut output_vector = FixedBitSet::with_capacity(vector.len());
        for bit_nr in 0..output_vector.len() {
            let mut enabled = true;
            // equal to: (vector & &self.exclusion_matrix[bit_nr]).is_clear()
            for (v, m) in zip_eq(vector.as_slice(), self.exclusion_matrix[bit_nr].as_slice()) {
                if (v & m) != 0u32 {
                    enabled = false;
                    break;
                }
            }
            output_vector.set(bit_nr, enabled);
        }
        output_vector
    }

    pub fn get_possibilities(&self, vector: &FixedBitSet) -> FixedBitSet {
        let non_singlet_mask: FixedBitSet = &self.tripoint_mask | &self.midpoint_mask;
        let possibilities: FixedBitSet = self.matrix_vector_multiply(vector);
        let masked_possibilities = &possibilities & &non_singlet_mask;

        if masked_possibilities.is_clear() {
            &possibilities & &self.singlet_mask
        } else {
            masked_possibilities
        }
    }

    pub fn get_bitarray(&self) -> FixedBitSet {
        self.filled_sites.clone()
    }

    pub fn solve(&self, find_all: bool) -> Vec<BitArraySolution> {
        let test_lattice = self.filled_sites.clone();

        let mut current_generation = HashSet::from([test_lattice]);
        let mut next_generation = HashSet::new();
        let mut depth = 0;
        let mut solutions = vec![];
        kdam::term::init(stderr().is_terminal());

        while solutions.is_empty() | (find_all & (next_generation.len() > 0)) {
            depth += 1;

            next_generation.clear();

            let iterator = tqdm!(
                current_generation.iter(),
                desc = format!("Current depth: {depth}"),
                mininterval = 1.0/60.0,
                bar_format = "{desc suffix=' '}|{animation}| {spinner} {count}/{total} [{percentage:.0}%] in {elapsed human=true} ({rate:.1}/s, eta: {remaining human=true})",
                colour = Colour::gradient(&["#5A56E0", "#EE6FF8"]),
                spinner = Spinner::new(
                    &["▁▂▃", "▂▃▄", "▃▄▅", "▄▅▆", "▅▆▇", "▆▇█", "▇█▇", "█▇▆", "▇▆▅", "▆▅▄", "▅▄▃", "▄▃▂", "▃▂▁", "▂▁▂"],
                    30.0,
                    1.0,
                )
            );
            for candidate in iterator {
                let possibilities = self.get_possibilities(candidate);

                if possibilities.is_clear() {
                    solutions.push(BitArraySolution(candidate.clone()));
                    continue;
                }

                for fillable_site in possibilities.ones() {
                    let mut new_candidate = candidate.clone();
                    new_candidate.set(fillable_site, true);

                    next_generation.insert(new_candidate);
                }
            }

            mem::swap(&mut current_generation, &mut next_generation);
        }
        solutions
    }

    pub fn __str__(&self) -> String {
        let mut output = String::from("BitArrayRepresentation: {\n");

        output += format!("  filled_sites = \n    {}\n", self.filled_sites).as_str();
        output += "  exclusion_matrix = \n";
        for (number, row) in self.exclusion_matrix.iter().enumerate() {
            output += format!("{number:5}. {}\n", row).as_str();
        }
        output += format!("  tripoint_mask = \n    {}\n", self.tripoint_mask).as_str();
        output += format!("  midpoint_mask = \n    {}\n", self.midpoint_mask).as_str();
        output += format!("  singlet_mask = \n    {}\n", self.singlet_mask).as_str();
        output += "}";

        output
    }

    pub fn __repr__(&self) -> String {
        format!("BitArrayRepresentation[{}]", self.filled_sites)
    }
}

pub struct Lattice {
    points: Vec<Arc<LatticePoint>>,
    oxygens: Vec<Oxygen>,
}

impl Lattice {
    /// Create an empty `Lattice`
    fn new() -> Self {
        Lattice {
            points: vec![],
            oxygens: vec![],
        }
    }

    /// Push a point to `self.points`
    fn add_point(&mut self, new_point: Arc<LatticePoint>) {
        self.points.push(new_point);
    }

    fn generate_exclusions(&mut self) {
        for (number, oxygen) in self.oxygens.iter().enumerate() {
            for index in oxygen.sitetype.iter() {
                let mut connections = self.points[index.0].get_connections().write().unwrap();
                connections.push(OxygenIndex(number));
            }
        }
        for oxygen in &mut self.oxygens {
            for connection in oxygen.sitetype.iter() {
                let point = &self.points[connection.0];
                for ox_con in point.get_connections().read().unwrap().iter() {
                    oxygen.exclusions.push(*ox_con);
                }
            }
            oxygen.exclusions.dedup();
        }
    }

    /// Turns the lattice problem in an abstracted form based on bitarrays.
    fn generate_intermediary(&self) -> BitArrayRepresentation {
        let filled_sites = FixedBitSet::with_capacity(self.oxygens.len());
        let mut exclusion_matrix: Vec<FixedBitSet> = vec![];

        let mut tripoint_mask = FixedBitSet::with_capacity(self.oxygens.len());
        let mut midpoint_mask = FixedBitSet::with_capacity(self.oxygens.len());
        let mut singlet_mask = FixedBitSet::with_capacity(self.oxygens.len());

        for (number, oxygen) in self.oxygens.iter().enumerate() {
            let mut exclusions = FixedBitSet::with_capacity(self.oxygens.len());
            for exclusion in &oxygen.exclusions {
                exclusions.set(exclusion.0, true);
            }
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
            tripoint_mask,
            midpoint_mask,
            singlet_mask,
        }
    }

    pub fn python_new(input_lattice: Vec<(Vec<f32>, Vec<Vec<f32>>)>) -> Self {
        let mut out_lattice = Lattice::new();
        // Create the silicon lattice
        for (location, ghosts) in input_lattice {
            let new_point = LatticePoint::new(location[0], location[1], None);
            out_lattice.add_point(new_point.clone());

            for ghost in ghosts {
                out_lattice.add_point(LatticePoint::new(
                    ghost[0],
                    ghost[1],
                    Some(new_point.clone()),
                ))
            }
        }
        out_lattice
            .points
            .sort_by_key(|p| (100.0 * p.x + p.y).round() as u32);

        // Fill in the oxygens
        let silicon_iterator = out_lattice.points.iter().map(|p| [p.x, p.y]).collect_vec();
        let kdtree: KdTree<_, 2> = (&silicon_iterator).into();
        let node_search_distance =
            kdtree.nearest_n::<SquaredEuclidean>(&[0f32, 0f32], 2)[1].distance * 1.1;

        // Tripoints
        for (number, silicon) in out_lattice.points.iter().enumerate() {
            let mut close_points =
                kdtree.within::<SquaredEuclidean>(&[silicon.x, silicon.y], node_search_distance);
            // Sort results on lattice number
            close_points.sort_by_key(|p| p.item);
            let sites = close_points
                .iter()
                .skip(1)
                .filter(|s| s.item as usize > number)
                .combinations(2)
                .filter(|a| {
                    out_lattice.points[a[0].item as usize].x
                        != out_lattice.points[a[1].item as usize].x
                });

            for site in sites {
                let x = (out_lattice.points[number].x
                    + out_lattice.points[site[0].item as usize].x
                    + out_lattice.points[site[1].item as usize].x)
                    / 3.0;
                let y = (out_lattice.points[number].y
                    + out_lattice.points[site[0].item as usize].y
                    + out_lattice.points[site[1].item as usize].y)
                    / 3.0;
                let sitetype = SiteType::Tripoint(Tripoint([
                    LatticeIndex(number),
                    LatticeIndex(site[0].item as usize),
                    LatticeIndex(site[1].item as usize),
                ]));
                out_lattice.oxygens.push(Oxygen::new(x, y, sitetype))
            }
        }

        // Midpoints
        for (number, silicon) in out_lattice.points.iter().enumerate() {
            let mut close_points =
                kdtree.within::<SquaredEuclidean>(&[silicon.x, silicon.y], node_search_distance);
            // Sort results on lattice number
            close_points.sort_by_key(|p| p.item);
            let sites = close_points
                .iter()
                .skip(1)
                .filter(|s| s.item as usize > number)
                .filter(|s| {
                    out_lattice.points[number].ghost_to.is_none()
                        || out_lattice.points[s.item as usize].ghost_to.is_none()
                });

            for site in sites {
                let x =
                    (out_lattice.points[number].x + out_lattice.points[site.item as usize].x) / 2.0;
                let y =
                    (out_lattice.points[number].y + out_lattice.points[site.item as usize].y) / 2.0;
                let sitetype = SiteType::Midpoint(Midpoint([
                    LatticeIndex(number),
                    LatticeIndex(site.item as usize),
                ]));
                out_lattice.oxygens.push(Oxygen::new(x, y, sitetype))
            }
        }

        // Singles
        for (number, silicon) in out_lattice
            .points
            .iter()
            .enumerate()
            .filter(|s| s.1.ghost_to.is_none())
        {
            out_lattice.oxygens.push(Oxygen::new(
                silicon.x,
                silicon.y,
                SiteType::Singlet(Singlet([LatticeIndex(number)])),
            ));
        }

        out_lattice.generate_exclusions();
        out_lattice
    }

    /// Returns the coordinates of the lattice points in two lists.
    /// Use with the * star operator in a `plt.plot` function:
    ///
    /// ```python
    /// plt.plot(*solved_lattice.points_to_plot(), "o")
    /// ```
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
    pub fn get_intermediary(&self) -> BitArrayRepresentation {
        self.generate_intermediary()
    }

    /// Returns a solved version of the lattice. Usefull for plotting
    /// and exporting.
    pub fn to_solved_lattice(&self, solution: &BitArraySolution) -> Self {
        let mut solved_oxygens = vec![];
        for number in 0..self.oxygens.len() {
            if solution.0[number] {
                solved_oxygens.push(self.oxygens[number].clone());
            }
        }
        Lattice {
            points: self.points.clone(),
            oxygens: solved_oxygens,
        }
    }

    /// Export the solved lattice to a json file.
    ///
    /// - path must be a valid path name.
    /// - name should end with ".json".
    pub fn export(&self, path: OsString, name: String) {
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

        let mut file = File::create(filename).unwrap();
        file.write_all(data.pretty(4).as_bytes()).unwrap();
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn new_mmx() {
//         let lattice_points = vec![
//             (vec![1.5, 0.0], vec![vec![1.5, 5.196152422706632]]),
//             (
//                 vec![0.0, 0.0],
//                 vec![
//                     vec![0.0, 5.196152422706632],
//                     vec![6.0, 0.0],
//                     vec![6.0, 5.196152422706632],
//                 ],
//             ),
//             (vec![2.25, 1.299038105676658], vec![]),
//             (
//                 vec![0.75, 1.299038105676658],
//                 vec![vec![6.75, 1.299038105676658]],
//             ),
//             (vec![1.5, 2.598076211353316], vec![]),
//             (
//                 vec![0.0, 2.598076211353316],
//                 vec![vec![6.0, 2.598076211353316]],
//             ),
//             (vec![2.25, 3.897114317029974], vec![]),
//             (
//                 vec![0.75, 3.897114317029974],
//                 vec![vec![6.75, 3.897114317029974]],
//             ),
//             (vec![4.5, 0.0], vec![vec![4.5, 5.196152422706632]]),
//             (vec![3.0, 0.0], vec![vec![3.0, 5.196152422706632]]),
//             (vec![5.25, 1.299038105676658], vec![]),
//             (vec![3.75, 1.299038105676658], vec![]),
//             (vec![4.5, 2.598076211353316], vec![]),
//             (vec![3.0, 2.598076211353316], vec![]),
//             (vec![5.25, 3.897114317029974], vec![]),
//             (vec![3.75, 3.897114317029974], vec![]),
//         ];
//         let lattice = Lattice::python_new(lattice_points);
//         let bit_lattice = lattice.get_intermediary();
//         let mut test_bitarray = bit_lattice.get_bitarray();
//         test_bitarray.set(0, true);

//         let old = dbg!(bit_lattice.get_possibilities_old(&test_bitarray));
//         println!("{old}");
//         let new = dbg!(bit_lattice.get_possibilities(&test_bitarray));
//         println!("{new}");
//         assert_eq!(old, new);
//     }
// }
