use std::process::exit;

use crystacean_rs::{bit_array_settings, BitArrayFilter, BitArraySettings, Lattice};

fn main() {
    let lattice = Lattice::from_dft_json("../test_lattices/T20.json".into(), 3.5, false);

    let options = bit_array_settings!(
        lattice,
        solve_filter = BitArrayFilter::Similarity,
        difference_distance = 0.1
    );
    // let filter = lattice.no_rings();
    let bit_lattice = lattice.get_intermediary(options);
    // let filtered = bit_lattice.filtered(filter);

    println!("{}", bit_lattice.__str__());
    bit_lattice.print_distances();
    exit(0);

    let solutions = bit_lattice.solve(true, false);
    println!(" Solutions found: {}", solutions.len());

    // println!("{}", filtered.__str__());

    // let filtered_solutions = filtered.solve(true, false);
    // println!(" Solutions found: {}", filtered_solutions.len());

    // let solutions_filtered = bit_lattice.solve_filtered(true, false);
    // println!("Solutions found: {}", solutions.len());
    // for (number, solution) in solutions.iter().enumerate() {
    //     let solution_lattice = lattice.to_solved_lattice(solution);
    //     solution_lattice.export(&OsString::from("../exports/test4"), format!("test2_{number}.json"));
    // }
}
