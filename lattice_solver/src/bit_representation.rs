use fixedbitset::FixedBitSet;

mod bit_rep_impl;
use crate::*;

pub struct BitArrayRepresentation {
    pub filled_sites: FixedBitSet,
    pub exclusion_matrix: Vec<FixedBitSet>,
    pub distances_matrix: Vec<Vec<f32>>,
    pub tripoint_mask: FixedBitSet,
    pub midpoint_mask: FixedBitSet,
    pub singlet_mask: FixedBitSet,
    pub filter: Option<FixedBitSet>,
    pub options: BitArraySettings,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BitArraySettings {
    pub max_singlets: usize,
    pub difference_distance: f32,
    pub max_x: f32,
    pub max_y: f32,
    pub solve_filter: BitArrayFilter,
}

impl BitArraySettings {
    pub const fn create(
        max_singlets: usize,
        difference_distance: f32,
        max: (f32, f32),
        solve_filter: BitArrayFilter,
    ) -> Self {
        Self {
            max_singlets,
            difference_distance,
            max_x: max.0,
            max_y: max.1,
            solve_filter,
        }
    }

    pub fn default(lattice: &Lattice) -> Self {
        let max = lattice.find_max();
        Self {
            max_singlets: 2,
            difference_distance: 0.05,
            max_x: max.0,
            max_y: max.1,
            solve_filter: BitArrayFilter::default(),
        }
    }
}

#[derive(Default)]
pub struct SettingsBuilder {
    max_singlets: Option<usize>,
    difference_distance: Option<f32>,
    max: Option<(f32, f32)>,
    solve_filter: Option<BitArrayFilter>,
}

impl SettingsBuilder {
    pub const fn max_singlets(self, value: usize) -> Self {
        Self {
            max_singlets: Some(value),
            difference_distance: self.difference_distance,
            max: self.max,
            solve_filter: self.solve_filter,
        }
    }
    pub const fn difference_distance(self, value: f32) -> Self {
        Self {
            max_singlets: self.max_singlets,
            difference_distance: Some(value),
            max: self.max,
            solve_filter: self.solve_filter,
        }
    }
    pub const fn max(self, value: (f32, f32)) -> Self {
        Self {
            max_singlets: self.max_singlets,
            difference_distance: self.difference_distance,
            max: Some(value),
            solve_filter: self.solve_filter,
        }
    }
    pub const fn solve_filter(self, value: BitArrayFilter) -> Self {
        Self {
            max_singlets: self.max_singlets,
            difference_distance: self.difference_distance,
            max: self.max,
            solve_filter: Some(value),
        }
    }
    pub fn build(self, lattice: &Lattice) -> BitArraySettings {
        BitArraySettings::create(
            self.max_singlets.unwrap_or(2),
            self.difference_distance.unwrap_or(0.05),
            self.max.unwrap_or_else(|| lattice.find_max()),
            self.solve_filter.unwrap_or_default(),
        )
    }
}

#[macro_export]
macro_rules! bit_array_settings {
    ( $latt:expr, $($setter_method: ident = $value: expr),*) => {
        // use lattice_solver::SettingsBuilder;
        lattice_solver::SettingsBuilder::default()$(.$setter_method($value))*.build(&$latt);
    };

    ( $latt:expr ) => {BitArraySettings::default(&$latt);};
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum BitArrayFilter {
    #[default]
    None,
    Similarity,
    SimTrees,
}

#[derive(Debug, PartialEq, Eq, Default)]
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

    #[must_use]
    pub fn __str__(&self) -> String {
        format!("{}", self.0)
    }
}
