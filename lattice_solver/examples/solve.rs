use crystacean_rs::{
    bit_array_settings, test_points::{huge_points, lattice_points}, BitArrayFilter, BitArraySettings, Lattice,
};

fn main() {
    let lattice = Lattice::python_new(huge_points(), 1.1, true);
    let tolerance = 0.9f32;

    let trees_options = bit_array_settings!(
        lattice,
        solve_filter = BitArrayFilter::Similarity,
        difference_distance = tolerance,
        max_singlets = 0
    );
    let trees_bit_lattice = lattice.get_intermediary(trees_options);

    // println!("{}", bit_lattice.__str__());

    let trees_solutions = trees_bit_lattice.solve(true, false);
    println!("\n\nSolutions found: {}", trees_solutions.len());

    // let old_options = bit_array_settings!(
    //     lattice,
    //     solve_filter = BitArrayFilter::Similarity,
    //     difference_distance = tolerance
    // );
    // let old_bit_lattice = lattice.get_intermediary(old_options);

    // // println!("{}", bit_lattice.__str__());

    // let old_solutions = old_bit_lattice.solve(true, false);
    // println!("\n\nSolutions found: {}", old_solutions.len());


    // let none_options = bit_array_settings!(
    //     lattice,
    //     solve_filter = BitArrayFilter::None,
    //     difference_distance = tolerance
    // );
    // let none_bit_lattice = lattice.get_intermediary(none_options);

    // // println!("{}", bit_lattice.__str__());

    // let none_solutions = none_bit_lattice.solve(true, false);
    // println!("\n\nSolutions found: {}", none_solutions.len());

    // dbg!(trees_solutions == old_solutions);
    // dbg!(trees_solutions == none_solutions);
    // let solutions_filtered = bit_lattice.solve_filtered(true, false);
    // println!("Solutions found: {}", solutions.len());
    // for (number, solution) in solutions.iter().enumerate() {
    //     let solution_lattice = lattice.to_solved_lattice(solution);
    //     solution_lattice.export(&OsString::from("../exports/test4"), format!("test2_{number}.json"));
    // }
}
