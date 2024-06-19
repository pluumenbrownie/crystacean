use std::collections::BTreeMap;

use lattice_solver::close_vector_tree::{CloseVectorTree, CloseVectorTreeMap};
use ordered_float::NotNan;

fn main() {
    let vector = vec![
        NotNan::new(1.1f32).unwrap(),
        NotNan::new(1.1f32).unwrap(),
        NotNan::new(1.5f32).unwrap(),
        NotNan::new(1.5f32).unwrap(),
        NotNan::new(2.5f32).unwrap(),
        NotNan::new(3.0f32).unwrap(),
    ];
    let tolerance = 0.4f32.try_into().unwrap();
    let mut map = dbg!(CloseVectorTreeMap::new(tolerance));

    dbg!(map.insert(vector, Default::default()));
    println!("{map:?}");

    let vector2 = vec![
        NotNan::new(1.1f32).unwrap(),
        NotNan::new(1.2f32).unwrap(),
        NotNan::new(1.4f32).unwrap(),
        NotNan::new(1.5f32).unwrap(),
        NotNan::new(2.3f32).unwrap(),
        NotNan::new(3.5f32).unwrap(),
    ];

    dbg!(map.insert(vector2, Default::default()));
    println!("{map:?}");
    
    let vector_fail = vec![
        NotNan::new(1.1f32).unwrap(),
        NotNan::new(1.2f32).unwrap(),
        NotNan::new(1.2f32).unwrap(),
        NotNan::new(1.5f32).unwrap(),
        NotNan::new(2.3f32).unwrap(),
        NotNan::new(3.5f32).unwrap(),
    ];

    dbg!(map.insert(vector_fail, Default::default()));
    println!("{map:?}");
    
    let vector_fail = vec![
        NotNan::new(1.1f32).unwrap(),
        NotNan::new(1.2f32).unwrap(),
        NotNan::new(1.2f32).unwrap(),
        NotNan::new(1.5f32).unwrap(),
        NotNan::new(2.0f32).unwrap(),
        NotNan::new(3.0f32).unwrap(),
    ];

    dbg!(map.insert(vector_fail, Default::default()));
    println!("{map:?}");
    
}
