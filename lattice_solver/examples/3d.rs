use std::{fs::File, io::Read};
use itertools::Itertools;
use json;


fn main() {
    let delta_y: f32 = 1.3;

    let mut buffer = String::new();
    let mut file = File::open("../DFT_results/T1S.json").unwrap();
    file.read_to_string(&mut buffer).unwrap();
    let parsed = json::parse(&buffer).unwrap();

    let last_id = &parsed["ids"][parsed["ids"].len() - 1].to_string();

    let numbers = &parsed[last_id]["positions"]["__ndarray__"][2];
    let atoms = numbers.members()
        .map(|j| j.as_f32().unwrap())
        .collect_vec()
        .chunks_exact(3)
        .map(|c| c.to_vec())
        .collect_vec();

    let ends = atoms.iter()
        .filter(|a| a[2] < 9.0)
        .collect_vec();

    for atom in ends {
        println!("{atom:?}");
    }
}