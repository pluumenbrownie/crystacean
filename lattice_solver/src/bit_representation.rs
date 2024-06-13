use fixedbitset::FixedBitSet;

mod bit_rep_impl;

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

#[derive(Debug, Clone, Copy)]
pub struct BitArraySettings {
    pub max_singlets: usize,
    pub difference_distance: f32,
    pub max_x: f32,
    pub max_y: f32,
}

impl BitArraySettings {
    pub const fn create(max_singlets: usize, difference_distance: f32, max: (f32, f32)) -> Self {
        Self {
            max_singlets,
            difference_distance,
            max_x: max.0,
            max_y: max.1,
        }
    }
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
