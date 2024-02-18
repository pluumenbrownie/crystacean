use std::{fs::File, io::Read};
use itertools::Itertools;
use lattice_solver::Lattice;


fn main() {
    // let delta_y: f32 = 1.4;

    // let mut buffer = String::new();
    // let mut file = File::open("../DFT_results/T1S.json").unwrap();
    // file.read_to_string(&mut buffer).unwrap();
    // let parsed = json::parse(&buffer).unwrap();

    // let last_id = &parsed["ids"][parsed["ids"].len() - 1].to_string();

    // let hydrogen_amount = &parsed[last_id]["numbers"]["__ndarray__"][2]
    //     .members()
    //     .filter_map(|n| n.as_usize())
    //     .filter(|n| n == &1)
    //     .count();
    // let numbers = &parsed[last_id]["positions"]["__ndarray__"][2];
    // let atoms = numbers.members()
    //     .map(|j| j.as_f32().unwrap())
    //     .collect_vec()
    //     .chunks_exact(3)
    //     .map(|c| c.to_vec())
    //     .collect_vec();

    // let ends = atoms.iter()
    //     .sorted_by(|&a, &b| a[2].total_cmp(&b[2]))
    //     .take(hydrogen_amount.div_ceil(2))
    //     .collect_vec();
    let lattice = Lattice::from_dft_json("../DFT_results/T1S.json".into(), 1.3);

    // for atom in ends {
    //     println!("{atom:?}");
    // }
}