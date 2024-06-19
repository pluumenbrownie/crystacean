use close_vector_tree::CloseVectorTreeMap;
use fixedbitset::FixedBitSet;
use itertools::{zip_eq, Itertools};
use kdam::{par_tqdm, tqdm, Colour, Spinner};
use ordered_float::NotNan;
use scc::Bag;
use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    io::{stderr, IsTerminal},
    mem,
};

use crate::*;
use rayon::prelude::*;

impl BitArrayRepresentation {
    /// Create a `BitArrayRepresentation` for testing purpouses.
    // #[cfg(doctest)]
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn create_debug(
        filled_sites: FixedBitSet,
        exclusion_matrix: Vec<FixedBitSet>,
        distances_matrix: Vec<Vec<f32>>,
        tripoint_mask: FixedBitSet,
        midpoint_mask: FixedBitSet,
        singlet_mask: FixedBitSet,
        filter: Option<FixedBitSet>,
        options: BitArraySettings,
    ) -> Self {
        Self {
            filled_sites,
            exclusion_matrix,
            distances_matrix,
            tripoint_mask,
            midpoint_mask,
            singlet_mask,
            filter,
            options,
        }
    }

    /// Performes a binary matrix-vector multiply between `self.exclusion_matrix`
    /// and the given solution `vector`. This method is used in `self.get_possibilities`
    /// to reveal the available silicon sites for a given solution.
    ///
    /// In this method, the input matrix is bitwise `AND`ed with each column of
    /// the exclusion matrix. The bits of the output vector are then set to 0
    /// if an AND operation between the input vector and the corresponding column
    /// of the matrix contain any ones.
    ///
    /// This method should not be used directly, use `self.get_possibilities` instead.
    ///
    /// ## Example:
    /// ```
    /// # use lattice_solver::{BitArraySolution, BitArrayRepresentation};
    /// # use fixedbitset::FixedBitSet;
    /// #
    /// let bit_array_repr = BitArrayRepresentation::create_debug(
    ///     //                    11000
    ///     //                    11000
    ///     // exclusion_matrix = 00110
    ///     //                    00110
    ///     //                    00001
    ///     # FixedBitSet::with_capacity_and_blocks(5, vec![0b00000]),
    ///     # vec![
    ///     #     FixedBitSet::with_capacity_and_blocks(5, vec![0b00011]),
    ///     #     FixedBitSet::with_capacity_and_blocks(5, vec![0b00011]),
    ///     #     FixedBitSet::with_capacity_and_blocks(5, vec![0b01100]),
    ///     #     FixedBitSet::with_capacity_and_blocks(5, vec![0b01100]),
    ///     #     FixedBitSet::with_capacity_and_blocks(5, vec![0b10000]),
    ///     # ],
    ///     # FixedBitSet::new(),
    ///     # FixedBitSet::new(),
    ///     # FixedBitSet::new(),
    ///     # None
    /// );
    ///
    /// let potential_solution = FixedBitSet::with_capacity_and_blocks(5, vec![0b11000]);
    /// let available_sites = bit_array_repr.matrix_vector_multiply(&potential_solution);
    ///
    /// let expected_answer = FixedBitSet::with_capacity_and_blocks(5, vec![0b00011]);
    /// assert_eq!(available_sites, expected_answer);
    /// ```
    /// Shown schematically:
    /// ```text
    /// 11000     0     1
    /// 11000     0     1
    /// 00110  x  0  =  0
    /// 00110     1     0
    /// 00001     1     0
    /// ```
    #[must_use]
    pub fn matrix_vector_multiply(&self, vector: &FixedBitSet) -> FixedBitSet {
        let mut output_vector = FixedBitSet::with_capacity(vector.len());
        for bit_nr in 0..output_vector.len() {
            output_vector.set(bit_nr, vector.is_disjoint(&self.exclusion_matrix[bit_nr]));
        }
        output_vector
    }

    /// Returns the valid available sites for the given `vector`.
    ///
    /// The valid sites are determined by:
    /// 1. Obtaining available sites by using `self.matrix_vector_multiply`.
    /// 1. available sites are checked with tripoint and midpoint mask.
    ///     - if tripoints or midpoints are available, singlet possibilities are
    ///       masked out before the results vector is returned.
    ///
    /// There are 5 output scenarios:
    /// # Outputs
    ///  - possibilties is empty                       -> `Ok(empty vector)`
    ///  - tri/mid `masked_possibilities` is empty
    ///      - available singlets <= `self.max_singlets` -> `Ok(vector)`
    ///  - tri/mid `masked_possibilities` is NOT empty   -> `Ok(masked_vector`)
    ///
    /// # Errors
    ///  - possibilities is empty after `rightmost_mask` -> `Err(())`
    ///  - tri/mid `masked_possibilities` is empty
    ///      - available singlets >  `self.max_singlets` -> `Err(())`
    /// Vectors which return `Err(())` should be ignored.
    ///
    /// ## Example:
    /// ```
    /// # use lattice_solver::{BitArraySolution, BitArrayRepresentation};
    /// # use fixedbitset::FixedBitSet;
    /// #
    /// let bit_array_repr = BitArrayRepresentation::create_debug(
    ///     //                    11000
    ///     //                    11000
    ///     // exclusion_matrix = 00110
    ///     //                    00110
    ///     //                    00001
    ///     # FixedBitSet::with_capacity_and_blocks(5, vec![0b00000]),
    ///     # vec![
    ///     #     FixedBitSet::with_capacity_and_blocks(5, vec![0b00011]),
    ///     #     FixedBitSet::with_capacity_and_blocks(5, vec![0b00011]),
    ///     #     FixedBitSet::with_capacity_and_blocks(5, vec![0b01100]),
    ///     #     FixedBitSet::with_capacity_and_blocks(5, vec![0b01100]),
    ///     #     FixedBitSet::with_capacity_and_blocks(5, vec![0b10000]),
    ///     # ],
    ///     // tripoint, midpoint and singlet masks:
    ///     FixedBitSet::with_capacity_and_blocks(5, vec![0b00011]),
    ///     FixedBitSet::with_capacity_and_blocks(5, vec![0b01100]),
    ///     FixedBitSet::with_capacity_and_blocks(5, vec![0b10000]),
    ///     # None
    /// );
    ///
    /// let tripoint_possible =
    ///     FixedBitSet::with_capacity_and_blocks(5, vec![0b00010]);
    /// let tripoint_possible_sites = bit_array_repr.get_possibilities(&tripoint_possible);
    /// let tripoint_possible_answer =
    ///     FixedBitSet::with_capacity_and_blocks(5, vec![0b01100]);
    /// assert_eq!(tripoint_possible_sites, Ok(tripoint_possible_answer));
    ///
    /// let tripoint_impossible =
    ///     FixedBitSet::with_capacity_and_blocks(5, vec![0b01001]);
    /// let tripoint_impossible_sites = bit_array_repr.get_possibilities(&tripoint_impossible);
    /// let tripoint_impossible_answer =
    ///     FixedBitSet::with_capacity_and_blocks(5, vec![0b10000]);
    /// assert_eq!(tripoint_impossible_sites, Ok(tripoint_impossible_answer));
    /// ```
    pub fn get_possibilities(&self, vector: &FixedBitSet) -> Result<FixedBitSet, &str> {
        let non_singlet_mask: FixedBitSet = &self.tripoint_mask | &self.midpoint_mask;
        let rightmost_bit = vector.maximum();

        // If no possibilities found here, vector is solution
        let mut possibilities: FixedBitSet = self.matrix_vector_multiply(vector);

        // possible_sites can be reused later
        let possible_sites = possibilities.count_ones(..);
        if possible_sites == 0 {
            return Ok(possibilities);
        };

        // Mask of bits covered in other threads. If a vector is empty after this, all
        // possible following states of this vector will be covered by other threads.
        if let Some(mask) = rightmost_bit {
            possibilities.set_range(..mask, false);
        };
        if possibilities.is_clear() {
            return Err("Vector is invalid.");
        };

        let masked_possibilities = &possibilities & &non_singlet_mask;
        if masked_possibilities.is_clear() {
            if possible_sites <= self.options.max_singlets {
                Ok(possibilities)
            } else {
                Err("Vector is invalid.")
            }
        } else {
            Ok(masked_possibilities)
        }
    }

    /// Returns a clone of `self.filled_sites`.
    #[must_use]
    pub fn get_bitarray(&self) -> FixedBitSet {
        self.filled_sites.clone()
    }

    #[must_use]
    pub fn filtered(&self, filter: site_filter::SiteFilter) -> Self {
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
        let mut distances_matrix = vec![];

        for (new_number, old_number) in filter_set.ones().enumerate() {
            tripoint_mask.set(new_number, self.tripoint_mask[old_number]);
            midpoint_mask.set(new_number, self.midpoint_mask[old_number]);
            singlet_mask.set(new_number, self.singlet_mask[old_number]);

            let mut new_matrix_row = FixedBitSet::with_capacity(new_length);
            let mut new_distances_row = vec![];
            for (col_number, old_col_number) in filter_set.ones().enumerate() {
                new_matrix_row.set(
                    col_number,
                    self.exclusion_matrix[old_number][old_col_number],
                );

                new_distances_row.push(self.distances_matrix[old_number][old_col_number]);
            }
            exclusion_matrix.push(new_matrix_row);
            distances_matrix.push(new_distances_row);
        }

        Self {
            filled_sites,
            exclusion_matrix,
            distances_matrix,
            tripoint_mask,
            midpoint_mask,
            singlet_mask,
            filter: Some(filter_set),
            options: self.options,
        }
    }

    // /// Starts the solving process.
    // ///
    // /// `find_all` can be set to `true` to find all
    // #[must_use]
    // pub fn solve(&self, find_all: bool, silent: bool) -> Vec<BitArraySolution> {
    //     let test_lattice = self.filled_sites.clone();

    //     let mut current_generation = vec![test_lattice];
    //     let mut next_generation = vec![];
    //     let mut depth = 0;
    //     let mut solutions = vec![];
    //     kdam::term::init(stderr().is_terminal());

    //     while solutions.is_empty() | (find_all & (!next_generation.is_empty())) {
    //         depth += 1;

    //         next_generation.clear();

    //         let iterator = if silent {
    //             tqdm!(
    //                 current_generation.iter(),
    //                 disable = true,
    //                 position = 1,
    //                 bar_format = ""
    //             )
    //         } else {
    //             tqdm!(
    //                 current_generation.iter(),
    //                 desc = format!("Current depth: {depth}"),
    //                 mininterval = 1.0/60.0,
    //                 bar_format = "{desc suffix=' '}|{animation}| {spinner} {count}/{total} [{percentage:.0}%] in {elapsed human=true} ({rate:.1}/s, eta: {remaining human=true})",
    //                 colour = Colour::gradient(&["#FF0000", "#FFDD00"]),
    //                 spinner = Spinner::new(
    //                     &["▁▂▃", "▂▃▄", "▃▄▅", "▄▅▆", "▅▆▇", "▆▇█", "▇█▇", "█▇▆", "▇▆▅", "▆▅▄", "▅▄▃", "▄▃▂", "▃▂▁", "▂▁▂"],
    //                     60.0,
    //                     1.0,
    //                 ),
    //                 leave = true,
    //                 position = 1
    //             )
    //         };

    //         for candidate in iterator {
    //             if let Ok(possibilities) = self.get_possibilities(candidate) {
    //                 if possibilities.is_clear() {
    //                     solutions.push(BitArraySolution(candidate.clone()));
    //                     continue;
    //                 }

    //                 for fillable_site in possibilities.ones() {
    //                     let mut new_candidate = candidate.clone();
    //                     new_candidate.set(fillable_site, true);

    //                     next_generation.push(new_candidate);
    //                 }
    //             }
    //         }

    //         mem::swap(&mut current_generation, &mut next_generation);
    //     }
    //     if let Some(filter) = &self.filter {
    //         for solution in &mut solutions {
    //             solution.inflate(filter);
    //         }
    //     }
    //     println!();
    //     solutions
    // }

    /// Starts the solving process.
    ///
    /// `find_all` can be set to `true` to find all
    #[must_use]
    pub fn solve(&self, find_all: bool, silent: bool) -> Vec<BitArraySolution> {
        let test_lattice = self.filled_sites.clone();

        let mut current_generation = vec![test_lattice];
        let mut next_generation = vec![];
        let mut depth = 0;
        let mut solutions = vec![];
        kdam::term::init(stderr().is_terminal());

        while solutions.is_empty() | (find_all & (!next_generation.is_empty())) {
            depth += 1;

            next_generation.clear();

            let iterator = if silent {
                tqdm!(
                    current_generation.iter(),
                    disable = true,
                    position = 1,
                    bar_format = ""
                )
            } else {
                tqdm!(
                    current_generation.iter(),
                    desc = format!("Current depth: {depth}"),
                    mininterval = 1.0/60.0,
                    bar_format = "{desc suffix=' '}|{animation}| {spinner} {count}/{total} [{percentage:.0}%] in {elapsed human=true} ({rate:.1}/s, eta: {remaining human=true})",
                    colour = Colour::gradient(&["#FF0000", "#FFDD00"]),
                    spinner = Spinner::new(
                        &["▁▂▃", "▂▃▄", "▃▄▅", "▄▅▆", "▅▆▇", "▆▇█", "▇█▇", "█▇▆", "▇▆▅", "▆▅▄", "▅▄▃", "▄▃▂", "▃▂▁", "▂▁▂"],
                        60.0,
                        1.0,
                    ),
                    leave = true,
                    position = 1
                )
            };

            let mut structure_map = HashMap::new();
            let mut new_structure_map =
                CloseVectorTreeMap::new(self.options.difference_distance.try_into().unwrap());

            for candidate in iterator {
                if let Ok(possibilities) = self.get_possibilities(candidate) {
                    if possibilities.is_clear() {
                        solutions.push(BitArraySolution(candidate.clone()));
                        continue;
                    }

                    for fillable_site in possibilities.ones() {
                        let mut new_candidate = candidate.clone();
                        new_candidate.set(fillable_site, true);

                        if self.solving_filter(
                            &new_candidate,
                            &mut structure_map,
                            &mut new_structure_map,
                        ) {
                            next_generation.push(new_candidate);
                        }
                    }
                }
            }

            mem::swap(&mut current_generation, &mut next_generation);
        }
        if let Some(filter) = &self.filter {
            for solution in &mut solutions {
                solution.inflate(filter);
            }
        }
        // println!();
        // println!();
        solutions
    }

    fn solving_filter(
        &self,
        new_candidate: &FixedBitSet,
        structure_map: &mut HashMap<(usize, usize, usize), Vec<Vec<f32>>>,
        new_structure_map: &mut CloseVectorTreeMap,
    ) -> bool {
        match self.options.solve_filter {
            BitArrayFilter::None => true,
            BitArrayFilter::Similarity => self.similarity_filter(new_candidate, structure_map),
            BitArrayFilter::New => self.new_similarity_filter(new_candidate, new_structure_map),
        }
    }

    pub(crate) fn similarity_filter(
        &self,
        new_candidate: &FixedBitSet,
        structure_map: &mut HashMap<(usize, usize, usize), Vec<Vec<f32>>>,
    ) -> bool {
        let type_distribution = (
            new_candidate.union_count(&self.tripoint_mask),
            new_candidate.union_count(&self.midpoint_mask),
            new_candidate.union_count(&self.singlet_mask),
        );

        let structure_vec: &Vec<Vec<f32>> = structure_map.entry(type_distribution).or_default();

        let mut new_structure = vec![];
        for (one, two) in new_candidate.ones().combinations(2).map(|v| (v[0], v[1])) {
            new_structure.push(self.distances_matrix[one][two]);
        }
        new_structure.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

        let unique = structure_vec.iter().all(|structure| {
            new_structure
                .iter()
                .zip(structure.iter())
                .any(|(one, two)| (one - two).abs() > self.options.difference_distance)
        });
        if unique {
            let structure_vec = structure_map.get_mut(&type_distribution).unwrap();
            structure_vec.push(new_structure);
        }
        unique
    }

    pub(crate) fn new_similarity_filter(
        &self,
        new_candidate: &FixedBitSet,
        structure_map: &mut CloseVectorTreeMap,
    ) -> bool {
        // let tolerance = self.options.difference_distance;
        // let type_distribution = (
        //     new_candidate.union_count(&self.tripoint_mask),
        //     new_candidate.union_count(&self.midpoint_mask),
        //     new_candidate.union_count(&self.singlet_mask),
        // );

        // let vec_of_trees: &Vec<BTreeMap<NotNan<f32>, Vec<usize>>> =
        //     structure_map.entry(type_distribution).or_default();

        // let new_structure = self.create_diff_vector(new_candidate);

        structure_map.insert(new_candidate, self)

        // let unique = if new_structure.len() == vec_of_trees.len() {
        //     new_structure
        //         .iter()
        //         .zip_eq(vec_of_trees)
        //         .map(|(value, tree)| {
        //             tree.range(value - tolerance..=value + tolerance)
        //                 .flat_map(|key| key.1)
        //                 .map(std::borrow::ToOwned::to_owned)
        //                 .collect()
        //         })
        //         .reduce(|acc: HashSet<usize>, h| acc.intersection(&h).copied().collect())
        //         .is_none()
        // } else {
        //     true
        // };

        // if unique {
        //     let vec_of_trees = structure_map.get_mut(&type_distribution).unwrap();
        //     if new_structure.len() == vec_of_trees.len() {
        //         for (value, tree) in zip_eq(new_structure.into_iter(), vec_of_trees.iter_mut()) {
        //             let new_id = tree.len();
        //             tree.entry(value).or_default().push(new_id);
        //         }
        //     } else {
        //         for value in new_structure {
        //             vec_of_trees.push(BTreeMap::new());
        //             vec_of_trees.last_mut().unwrap().insert(value, vec![0]);
        //         }
        //     }
        // }
        // unique
    }

    pub fn create_diff_vector(&self, new_candidate: &FixedBitSet) -> Vec<NotNan<f32>> {
        let mut new_structure = vec![];
        for (one, two) in new_candidate.ones().combinations(2).map(|v| (v[0], v[1])) {
            new_structure.push(NotNan::new(self.distances_matrix[one][two]).unwrap());
        }
        new_structure.sort_unstable();
        new_structure
    }

    /// Starts the solving process.
    ///
    /// `find_all` can be set to `true` to find all
    ///
    /// # Panics
    /// Could technically panic but I don't see that happening.
    #[must_use]
    pub fn solve_parallel(&self, find_all: bool, silent: bool) -> Vec<BitArraySolution> {
        let test_lattice = self.filled_sites.clone();

        let mut current_generation = vec![test_lattice];
        let mut next_generation = vec![];
        let mut depth = 0;
        // let mut solutions = Mutex::new(vec![]);
        let solution_bag: Bag<BitArraySolution> = Bag::default();
        kdam::term::init(stderr().is_terminal());

        while solution_bag.is_empty() | (find_all & (!next_generation.is_empty())) {
            depth += 1;

            next_generation.clear();

            let iterator = if silent {
                par_tqdm!(
                    current_generation.par_iter(),
                    disable = true,
                    position = 1,
                    bar_format = ""
                )
            } else {
                par_tqdm!(
                    current_generation.par_iter(),
                    desc = format!("Current depth: {depth}"),
                    mininterval = 1.0/60.0,
                    bar_format = "{desc suffix=' '}|{animation}| {spinner} {count}/{total} [{percentage:.0}%] in {elapsed human=true} ({rate:.1}/s, eta: {remaining human=true})",
                    colour = Colour::gradient(&["#0000FF", "#00FFFF"]),
                    spinner = Spinner::new(
                        &["▁▂▃", "▂▃▄", "▃▄▅", "▄▅▆", "▅▆▇", "▆▇█", "▇█▇", "█▇▆", "▇▆▅", "▆▅▄", "▅▄▃", "▄▃▂", "▃▂▁", "▂▁▂"],
                        60.0,
                        1.0,
                    ),
                    leave = true,
                    position = 1
                )
            };

            // let structure_archive: scc::HashMap<(usize, usize, usize), Vec<Vec<f32>>> =
            //     scc::HashMap::new();

            next_generation = iterator
                .map(|vector| (vector, self.get_possibilities(vector)))
                .filter_map(|(v, c)| c.map_or_else(|_| None, |p| Some((v, p))))
                .map(|(v, c)| {
                    let mut new_candidates = vec![];
                    for fillable_site in c.ones() {
                        let mut new = v.clone();
                        new.set(fillable_site, true);
                        new_candidates.push(new);
                    }
                    (v, new_candidates)
                })
                .flat_map_iter(|(v, c)| {
                    if c.is_empty() {
                        solution_bag.push(BitArraySolution(v.clone()));
                    }
                    c
                })
                .collect();

            mem::swap(&mut current_generation, &mut next_generation);
        }
        let mut solutions = vec![];
        for mut solution in solution_bag {
            if let Some(filter) = &self.filter {
                solution.inflate(filter);
            }
            solutions.push(solution);
        }
        solutions
    }

    pub fn print_distances(&self) {
        for row in &self.distances_matrix {
            println!("{row:?}");
        }
    }

    #[must_use]
    pub fn __str__(&self) -> String {
        let mut output = String::from("BitArrayRepresentation: {\n");

        output += format!("  filled_sites = \n    {}\n", self.filled_sites).as_str();
        output += "  exclusion_matrix = \n";
        for (number, row) in self.exclusion_matrix.iter().enumerate() {
            output += format!("{number:5}. {row}\n").as_str();
        }
        output += format!("  tripoint_mask = \n    {}\n", self.tripoint_mask).as_str();
        output += format!("  midpoint_mask = \n    {}\n", self.midpoint_mask).as_str();
        output += format!("  singlet_mask = \n    {}\n", self.singlet_mask).as_str();
        output += "  filter = \n    ";
        output += self
            .filter
            .as_ref()
            .map_or_else(|| "None".into(), |fbs| format!("{fbs}"))
            .as_str();
        output += "\n}";

        output
    }

    #[must_use]
    pub fn __repr__(&self) -> String {
        format!("BitArrayRepresentation[{}]", self.filled_sites)
    }
}
