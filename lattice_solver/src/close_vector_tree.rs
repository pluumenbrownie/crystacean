use std::collections::{BTreeMap, HashMap, HashSet};

use fixedbitset::FixedBitSet;
use fmt_derive;
use itertools::Itertools;
use ordered_float::NotNan;

use crate::BitArrayRepresentation;

#[derive(Debug)]
pub struct CloseVectorTree {
    size: usize,
    tolerance: NotNan<f32>,
    vector: Vec<BTreeMap<NotNan<f32>, Vec<usize>>>,
}

impl CloseVectorTree {
    fn pop_size(&mut self) -> usize {
        let answer = self.size;
        self.size += 1;
        answer
    }

    pub fn length(length: usize, tolerance: NotNan<f32>) -> Self {
        let mut vector = vec![];
        for _ in 0..length {
            vector.push(BTreeMap::default());
        }
        Self {
            vector,
            tolerance,
            size: 0,
        }
    }

    pub fn check(&self, vector: &Vec<NotNan<f32>>) -> bool {
        let mut per_value = vec![];
        for (value, tree) in vector.iter().zip_eq(&self.vector) {
            let in_range: Vec<usize> = tree
                .range(value - self.tolerance..=value + self.tolerance)
                .flat_map(|(_, v)| v)
                .dedup()
                .copied()
                .collect();
            // println!("{in_range:?}");
            per_value.push(in_range);
        }
        let alike_vectors = per_value
            .into_iter()
            // .reduce(|acc, h| acc.intersection(&h).copied().collect())
            .reduce(|acc, h| {
                acc.into_iter()
                    .filter(|a| h.contains(a))
                    .collect_vec()
            })
            .unwrap_or_default();
        // println!("{alike_vectors:?}");
        alike_vectors.is_empty()
    }

    fn insert_blind(&mut self, vector: Vec<NotNan<f32>>) {
        let id = self.pop_size();
        for (&key, tree) in vector.iter().zip_eq(&mut self.vector) {
            tree.entry(key).or_insert(vec![]).push(id);
        }
    }

    pub fn insert(&mut self, vector: Vec<NotNan<f32>>) -> bool {
        if self.check(&vector) {
            self.insert_blind(vector);
            true
        } else {
            false
        }
    }
}

#[derive(fmt_derive::Debug)]
// pub struct CloseVectorTreeMap {
pub struct CloseVectorTreeMap {
    map: HashMap<PointTypeCount, CloseVectorTree>,
    tolerance: NotNan<f32>,
}

impl CloseVectorTreeMap {
    // pub fn insert(&mut self, vector: Vec<NotNan<f32>>) -> bool {
    //     let mut close_vector_tree = self.map.entry();
    //     todo!()
    // }
    pub fn new(tolerance: NotNan<f32>) -> Self {
        Self {
            map: HashMap::default(),
            tolerance,
        }
    }

    pub fn get(&mut self, point_type_count: PointTypeCount, length: usize) -> &CloseVectorTree {
        self.map
            .entry(point_type_count)
            .or_insert_with(|| CloseVectorTree::length(length, self.tolerance))
    }

    pub fn get_mut(
        &mut self,
        point_type_count: PointTypeCount,
        length: usize,
    ) -> &mut CloseVectorTree {
        self.map
            .entry(point_type_count)
            .or_insert_with(|| CloseVectorTree::length(length, self.tolerance))
    }

    pub fn insert(&mut self, bitvector: &FixedBitSet, repr: &BitArrayRepresentation) -> bool {
        let point_type_count = PointTypeCount::create(bitvector, repr);
        let vector = repr.create_diff_vector(bitvector);
        let tree = self.get_mut(point_type_count, vector.len());
        tree.insert(vector)
    }
}

#[derive(Hash, PartialEq, Eq, Debug, Default)]
pub struct PointTypeCount {
    tri: usize,
    mid: usize,
    sin: usize,
}

impl PointTypeCount {
    fn create(bitvector: &FixedBitSet, repr: &BitArrayRepresentation) -> Self {
        Self {
            tri: bitvector.union_count(&repr.tripoint_mask),
            mid: bitvector.union_count(&repr.midpoint_mask),
            sin: bitvector.union_count(&repr.singlet_mask),
        }
    }
}
