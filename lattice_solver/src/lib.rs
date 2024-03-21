// Optimisation Ideas:
// - Turn current/next generation into hashsets to remove archive requirement
//     - wont work with depth first
// - Try depth first and speed that up
// - parallellise
//     - try rayon
// - get profiler to work

use fixedbitset::FixedBitSet;
use itertools::{izip, zip_eq, Itertools};
use json::{object, JsonValue};
use kdam::{tqdm, Colour, Spinner};
use std::{
    collections::HashSet, f32::consts::PI, ffi::OsString, fs::File, io::{stderr, IsTerminal}, iter::zip, mem, sync::{Arc, RwLock}
};

use kiddo::{KdTree, SquaredEuclidean};
use std::io::prelude::*;

const TRIPLE_CROWN: [f32; 9] = [
    -0.000082766,
    -1.397718937,
    -0.517599594,

    -1.222253196,
    0.702780512,
    -0.517599594,

    1.204777673,
    0.711865236,
    -0.517599594,
];
const DOUBLE_CROWN: [f32; 6] = [
    2.2252961107321143 - 1.0526814404964109,
    -1.3867293215799141 - -1.0917258490752173,
    -(9.26392293591872 - 8.41953003801284),

    2.2252961107321143 - 3.3223088674933625,
    -1.3867293215799141 - -1.6248356621955473,
    -(9.26392293591872 - 8.313286837087986),
];
const SINGLE_CROWN: [f32; 3] = [0.0, 0.0, -1.7];

fn double_crown_rotated(theta: f32) -> [f32; 6] {
    // println!("{theta:?}");
    let r = ((2.2252961107321143 - 1.0526814404964109 as f32).powi(2) + (-1.3867293215799141 - -1.0917258490752173 as f32).powi(2)).sqrt();
    [
        r*(theta).cos(),
        r*(theta).sin(),
        -(9.26392293591872 - 8.41953003801284),

        r*(theta + PI).cos(),
        r*(theta + PI).sin(),
        -(9.26392293591872 - 8.313286837087986),
    ]
}

fn triple_crown_rotated(theta: f32) -> [f32; 9] {
    // println!("{theta:?}");
    let r = ((-0.000082766 as f32).powi(2) + (-1.397718937 as f32).powi(2)).sqrt();
    [
        r*(theta).cos(),
        r*(theta).sin(),
        -0.517599594,

        r*(theta + (2.0/3.0 * PI)).cos(),
        r*(theta + (2.0/3.0 * PI)).sin(),
        -0.517599594,

        r*(theta - (2.0/3.0 * PI)).cos(),
        r*(theta - (2.0/3.0 * PI)).sin(),
        -0.517599594,
    ]
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct OxygenIndex(usize);

#[derive(Clone, Copy, Debug)]
struct LatticeIndex(usize);

#[derive(Debug)]
struct LatticePoint {
    x: f32,
    y: f32,
    z: f32,
    connected_to: RwLock<Vec<OxygenIndex>>,
    ghost_to: Option<Arc<LatticePoint>>,
}

impl LatticePoint {
    fn new(x: f32, y: f32, z: f32, ghost_to: Option<Arc<LatticePoint>>) -> Arc<Self> {
        Arc::new(LatticePoint {
            x,
            y,
            z,
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

    fn distance_squared_to(&self, other: &LatticePoint) -> f32 {
        (self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2)
    }
}

#[derive(Clone)]
struct Oxygen {
    x: f32,
    y: f32,
    z: f32,
    sitetype: SiteType,
    exclusions: Vec<OxygenIndex>,
}

impl Oxygen {
    fn new(x: f32, y: f32, z: f32, sitetype: SiteType) -> Self {
        Oxygen {
            x,
            y,
            z,
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

#[derive(Debug, PartialEq)]
pub struct BitArraySolution(pub FixedBitSet);

impl BitArraySolution {
    /// Convert a solution from a filtered `BitArrayRepresentation` back into
    /// it's full form.
    ///
    /// ```
    /// use lattice_solver::BitArraySolution;
    /// use fixedbitset::FixedBitSet;
    ///
    /// let mut compressed = BitArraySolution(
    ///     FixedBitSet::with_capacity_and_blocks(4, vec![0b10110])
    /// );
    /// let full = BitArraySolution(
    ///     FixedBitSet::with_capacity_and_blocks(6, vec![0b1001010000])
    /// );
    /// let filter = FixedBitSet::with_capacity_and_blocks(6, vec![0b1101011000]);
    ///
    /// compressed.inflate(&filter);
    /// assert_eq!(compressed, full);
    /// ```
    /// Shown schematically:
    /// ```text
    /// 10 1 10
    /// 1101011000
    /// 1001010000
    /// ```
    pub fn inflate(&mut self, filter_bitset: &FixedBitSet) {
        let mut new_vector = FixedBitSet::with_capacity(filter_bitset.len());
        for (self_number, new_location) in filter_bitset.ones().enumerate() {
            new_vector.set(new_location, self.0[self_number]);
        }
        self.0 = new_vector;
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.0)
    }
}

pub struct BitArrayRepresentation {
    filled_sites: FixedBitSet,
    exclusion_matrix: Vec<FixedBitSet>,
    tripoint_mask: FixedBitSet,
    midpoint_mask: FixedBitSet,
    singlet_mask: FixedBitSet,
    filter: Option<FixedBitSet>,
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

    pub fn filtered(&self, filter: SiteFilter) -> BitArrayRepresentation {
        // println!("{filter:?}");
        let mut filter_set = FixedBitSet::with_capacity(self.filled_sites.len());
        filter_set.toggle_range(..);
        for number in filter.wrapped {
            filter_set.set(number.0, false);
        }
        let new_length = filter_set.count_ones(..);

        let filled_sites = FixedBitSet::with_capacity(new_length);

        let mut tripoint_mask = FixedBitSet::with_capacity(new_length);
        let mut midpoint_mask = FixedBitSet::with_capacity(new_length);
        let mut singlet_mask = FixedBitSet::with_capacity(new_length);

        let mut exclusion_matrix = vec![];

        for (new_number, old_number) in filter_set.ones().enumerate() {
            tripoint_mask.set(new_number, self.tripoint_mask[old_number]);
            midpoint_mask.set(new_number, self.midpoint_mask[old_number]);
            singlet_mask.set(new_number, self.singlet_mask[old_number]);

            let mut new_matrix_row = FixedBitSet::with_capacity(new_length);
            for (col_number, old_col_number) in filter_set.ones().enumerate() {
                new_matrix_row.set(
                    col_number,
                    self.exclusion_matrix[old_number][old_col_number].into(),
                );
            }
            exclusion_matrix.push(new_matrix_row);
        }

        BitArrayRepresentation {
            filled_sites,
            exclusion_matrix,
            tripoint_mask,
            midpoint_mask,
            singlet_mask,
            filter: Some(filter_set),
        }
    }

    pub fn solve(&self, find_all: bool) -> Vec<BitArraySolution> {
        let test_lattice = self.filled_sites.clone();

        let mut current_generation = HashSet::from([test_lattice]);
        let mut next_generation = HashSet::new();
        let mut depth = 0;
        let mut solutions = vec![];
        kdam::term::init(stderr().is_terminal());

        while solutions.is_empty() | (find_all & (!next_generation.is_empty())) {
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
                    60.0,
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
        if let Some(filter) = &self.filter {
            for solution in &mut solutions {
                solution.inflate(filter);
            }
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
        output += format!("  filter = \n    ").as_str();
        output += match &self.filter {
            None => "None".into(),
            Some(fbs) => format!("{}", fbs),
        }
        .as_str();
        output += "\n}";

        output
    }

    pub fn __repr__(&self) -> String {
        format!("BitArrayRepresentation[{}]", self.filled_sites)
    }
}

#[derive(Clone, Debug)]
pub struct SiteFilter {
    wrapped: Vec<OxygenIndex>,
}

impl SiteFilter {
    pub fn empty() -> Self {
        SiteFilter { wrapped: vec![] }
    }
}

pub struct Lattice {
    points: Vec<Arc<LatticePoint>>,
    oxygens: Vec<Oxygen>,
    source_file: Option<JsonValue>,
}

impl Lattice {
    /// Create an empty `Lattice`
    fn new() -> Self {
        Lattice {
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
            filter: None,
        }
    }

    /// distance_margin should be 1.1 for 2D, 1.4 for 3D
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
        let kdtree: KdTree<_, 3> = (&silicon_iterator).into();
        let node_search_distance = if autodetect_margin {
            kdtree.nearest_n::<SquaredEuclidean>(&first_point_location, 2)[1].distance
                * distance_margin
        } else {
            distance_margin.powi(2)
        };
        // println!("{node_search_distance:?}");

        // Tripoints
        let mut covered_sites = HashSet::new();
        for (number, silicon) in out_lattice.points.iter().enumerate() {
            let mut close_points = kdtree.within::<SquaredEuclidean>(
                &[silicon.x, silicon.y, silicon.z],
                node_search_distance,
            );
            // Sort results on lattice number
            close_points.sort_by_key(|p| p.item);
            let sites = close_points
                .iter()
                .filter(|s| s.item as usize != number)
                .combinations(2)
                .filter(|a| {
                    let mut identifier = [number as u64, a[0].item, a[1].item];
                    identifier.sort();
                    covered_sites.insert(identifier)
                })
                .filter(|a| {
                    out_lattice.points[a[0].item as usize].x
                        != out_lattice.points[a[1].item as usize].x
                })
                .filter(|a| {
                    out_lattice.points[number].ghost_to.is_none()
                        || out_lattice.points[a[0].item as usize].ghost_to.is_none()
                        || out_lattice.points[a[1].item as usize].ghost_to.is_none()
                })
                .filter(|a| {
                    out_lattice.points[a[0].item as usize]
                        .distance_squared_to(&out_lattice.points[a[1].item as usize])
                        <= node_search_distance
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
                let z = (out_lattice.points[number].z
                    + out_lattice.points[site[0].item as usize].z
                    + out_lattice.points[site[1].item as usize].z)
                    / 3.0
                    - 1.1;
                let sitetype = SiteType::Tripoint(Tripoint([
                    LatticeIndex(number),
                    LatticeIndex(site[0].item as usize),
                    LatticeIndex(site[1].item as usize),
                ]));
                out_lattice.oxygens.push(Oxygen::new(x, y, z, sitetype))
            }
        }

        // Midpoints
        for (number, silicon) in out_lattice.points.iter().enumerate() {
            let mut close_points = kdtree.within::<SquaredEuclidean>(
                &[silicon.x, silicon.y, silicon.z],
                node_search_distance,
            );
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
                let z = (out_lattice.points[number].z + out_lattice.points[site.item as usize].z)
                    / 2.0
                    - 1.4;
                let sitetype = SiteType::Midpoint(Midpoint([
                    LatticeIndex(number),
                    LatticeIndex(site.item as usize),
                ]));
                out_lattice.oxygens.push(Oxygen::new(x, y, z, sitetype))
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
                silicon.z - 1.7,
                SiteType::Singlet(Singlet([LatticeIndex(number)])),
            ));
        }

        out_lattice.generate_exclusions();
        out_lattice
    }

    fn add_source_file(&mut self, source_file: JsonValue) {
        self.source_file = Some(source_file);
    }

    pub fn from_dft_json(filename: String, distance_margin: f32, autodetect_margin: bool) -> Self {
        let mut buffer = String::new();
        let mut file = File::open(filename).expect("Opening file failed.");
        file.read_to_string(&mut buffer)
            .expect("Reading file failed.");
        let parsed = json::parse(&buffer).expect("Parsing file failed/");

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
            .map(|c| c.to_vec())
            .zip(parsed[last_id]["numbers"]["__ndarray__"][2].members().map(|i| i.as_i32()))
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
            .expect("Json 'cell' property format is incorrect.");

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

        let mut lattice = Lattice::python_new(input_lattice, distance_margin, autodetect_margin);
        lattice.add_source_file(parsed.clone());

        lattice
    }

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

    pub fn no_rings_plot(&self) -> Vec<(usize, f32, f32)> {
        let parsed = self.source_file.as_ref().unwrap();

        let last_id = &parsed["ids"][parsed["ids"].len() - 1].to_string();

        let numbers = &parsed[last_id]["positions"]["__ndarray__"][2];
        let atoms = numbers
            .members()
            .map(|j| j.as_f32().unwrap())
            .collect_vec()
            .chunks_exact(3)
            .map(|c| c.to_vec())
            .collect_vec();

        let top_silicon_locations = atoms
            .iter()
            .sorted_by(|&a, &b| a[2].total_cmp(&b[2]))
            .filter(|&s| (s[2] > 8.4) & (s[2] < 10.0))
            .map(|v| [v[0], v[1], v[2]])
            .collect_vec();

        let points_vector = self.points.iter().map(|p| [p.x, p.y, p.z]).collect_vec();
        println!("{points_vector:?}");
        println!("length: {}", points_vector.len());
        let points_tree: KdTree<_, 3> = (&points_vector).into();

        let point_group_vector = {
            let mut point_group_vector = vec![1000; self.points.len()];

            for (number, atom) in top_silicon_locations.iter().enumerate() {
                let close_points = points_tree.within::<SquaredEuclidean>(atom, 1.7f32.powi(2));
                println!("{:?} --- {:?}", atom, close_points);
                for point in close_points {
                    point_group_vector[point.item as usize] = number;
                }
            }
            point_group_vector
        };
        println!("{point_group_vector:?}");

        izip!(
            point_group_vector,
            self.points.iter().map(|o| o.x),
            self.points.iter().map(|o| o.y)
        )
        .collect_vec()
    }

    pub fn no_rings(&self) -> SiteFilter {
        let parsed = self.source_file.as_ref().unwrap();

        let last_id = &parsed["ids"][parsed["ids"].len() - 1].to_string();

        let numbers = &parsed[last_id]["positions"]["__ndarray__"][2];
        let atoms = numbers
            .members()
            .map(|j| j.as_f32().unwrap())
            .collect_vec()
            .chunks_exact(3)
            .map(|c| c.to_vec())
            .collect_vec();

        let top_silicon_locations = atoms
            .iter()
            .sorted_by(|&a, &b| a[2].total_cmp(&b[2]))
            .filter(|&s| (s[2] > 8.4) & (s[2] < 10.0))
            .map(|v| [v[0], v[1], v[2]])
            .collect_vec();

        let points_vector = self.points.iter().map(|p| [p.x, p.y, p.z]).collect_vec();
        // println!("{points_vector:?}");
        // println!("length: {}", points_vector.len());
        let points_tree: KdTree<_, 3> = (&points_vector).into();

        let point_group_vector = {
            let mut point_group_vector = vec![1000; self.points.len()];

            for (number, atom) in top_silicon_locations.iter().enumerate() {
                let close_points = points_tree.within::<SquaredEuclidean>(atom, 1.7f32.powi(2));
                // println!("{:?} --- {:?}", atom, close_points);
                for point in close_points {
                    point_group_vector[point.item as usize] = number;
                }
            }
            point_group_vector
        };
        // println!("{point_group_vector:?}");

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

        SiteFilter {
            wrapped: disabled_oxygens,
        }
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
            source_file: self.source_file.clone(),
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

        let mut file = File::create(&filename).expect(&*format!("{:?}", &filename));
        file.write_all(data.pretty(4).as_bytes()).unwrap();
    }

    pub fn export_as_ase_json(&self, filename: String) {
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
        };

        let mut file = File::create(&filename).expect("Folder does not exist!");
        file.write_all(export_data.pretty(4).as_bytes()).unwrap();
        println!("Saved {filename}");
    }

    fn add_crown(&self, oxygen: &Oxygen, new_numbers: &mut JsonValue, new_positions: &mut JsonValue) {
        let double_crown;
        let single_crown;
        let new_crown = match oxygen.sitetype {
            SiteType::Singlet(_) => {
                single_crown = triple_crown_rotated(0.0);
                single_crown.iter()
            },
            SiteType::Midpoint(p) => {
                double_crown = double_crown_rotated(
                    double_angle(&self.points[p.0[0].0], &self.points[p.0[1].0])
                );
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

fn double_angle(p1: &LatticePoint, p2: &LatticePoint) -> f32 {
    PI - ((p2.y-p1.y) / ((p2.x-p1.x).powi(2)+(p2.y-p1.y).powi(2)).sqrt() ).acos()
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
            ))
        }
    }
    out_lattice
        .points
        .sort_by_key(|p| (100.0 * p.x + p.y).round() as u32);
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
