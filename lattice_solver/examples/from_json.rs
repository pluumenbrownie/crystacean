#![allow(clippy::excessive_precision)]

use lattice_solver::{bit_array_settings, BitArrayFilter, BitArraySettings, Lattice};

fn main() {
    let lattice = Lattice::from_dft_json("../exports/16_base_2.json".into(), 1.1, true);

    let options = bit_array_settings!(lattice);
    let bit_lattice = lattice.get_intermediary(options);

    // bit_lattice.print_distances();
    println!("{}", bit_lattice.__str__());

    let solutions = bit_lattice.solve(true, false);
    println!("Solutions found: {}", solutions.len());

    // let solutions_filtered = bit_lattice.solve_filtered(true, false);
    // println!("Solutions found: {}", solutions.len());
    // for (number, solution) in solutions.iter().enumerate() {
    //     let solution_lattice = lattice.to_solved_lattice(solution);
    //     solution_lattice.export(&OsString::from("../exports/test4"), format!("test2_{number}.json"));
    // }
}
