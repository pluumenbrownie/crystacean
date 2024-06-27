use lattice_solver::{bit_array_settings, BitArraySettings, Lattice};

fn main() {
    let lattice = Lattice::from_dft_json("../exports/json_tester.json".into(), 5.0, false);

    let options = bit_array_settings!(lattice, max_singlets = 0);
    let filter = lattice.no_rings();
    let bit_lattice = lattice.get_intermediary(options);
    let filtered = bit_lattice.filtered(filter);

    // bit_lattice.print_distances();
    println!("{}", bit_lattice.__str__());

    let solutions = bit_lattice.solve(true, false);
    println!(" Solutions found: {}", solutions.len());

    println!("{}", filtered.__str__());

    let filtered_solutions = filtered.solve(true, false);
    println!(" Solutions found: {}", filtered_solutions.len());

    // let solutions_filtered = bit_lattice.solve_filtered(true, false);
    // println!("Solutions found: {}", solutions.len());
    // for (number, solution) in solutions.iter().enumerate() {
    //     let solution_lattice = lattice.to_solved_lattice(solution);
    //     solution_lattice.export(&OsString::from("../exports/test4"), format!("test2_{number}.json"));
    // }
}
