use lattice_solver::{
    bit_array_settings,
    test_points::{huge_points, lattice_points},
    BitArrayFilter, BitArraySettings, BitArraySolution, Lattice,
};
use termion::{clear, cursor};

// #[divan::bench(max_time = 120, args = [true, false])]
// fn full_run_parallel(silent: bool) -> Vec<lattice_solver::BitArraySolution> {
//     let lattice_points = lattice_points();

//     let lattice = Lattice::python_new(lattice_points, 1.1, true);

//     let options = bit_array_settings!(lattice);
//     let bit_lattice = lattice.get_intermediary(options);

//     let result = bit_lattice.solve_parallel(true, silent);
//     print!("\r{}{}", cursor::Right(22), clear::AfterCursor);
//     result
// }

#[divan::bench(max_time = 120, args = [BitArrayFilter::Similarity, BitArrayFilter::Flipped, BitArrayFilter::InsideOut])]
#[ignore]
fn huge_run(filter: BitArrayFilter) -> Vec<BitArraySolution> {
    println!();
    let lattice = Lattice::from_dft_json(
        "../exports/scaling_lattices/T16.json".into(),
        1.1,
        true,
    );

    let options = bit_array_settings!(
        lattice,
        solve_filter = filter,
        difference_distance = 0.05,
        max_singlets = 0
    );

    let bit_lattice = lattice.get_intermediary(options);

    let result = bit_lattice.solve(true, false);
    print!(
        "\r{}{}{}",
        cursor::Up(1),
        cursor::Right(19),
        clear::AfterCursor
    );
    result
}

#[divan::bench(max_time = 120, args = ["T04.json", "T08.json", "T09.json", "T12.json", "T12alt.json", "T16.json", "T20.json"])]
// #[ignore]
fn scaling(structure: &str) -> Vec<BitArraySolution> {
    let lattice = Lattice::from_dft_json(
        format!("../exports/scaling_lattices/{}", structure),
        1.1,
        true,
    );

    let options = bit_array_settings!(
        lattice,
        solve_filter = BitArrayFilter::None,
        difference_distance = 0.1,
        max_singlets = 0
    );

    let bit_lattice = lattice.get_intermediary(options);

    let result = bit_lattice.solve(true, true);
    print!(
        "\r{}{}",
        cursor::Right(19),
        clear::AfterCursor
    );
    result
}

#[divan::bench(max_time = 120, args = ["T04.json", "T08.json", "T09.json", "T12.json", "T12alt.json", "T16.json", "T20.json"])]
// #[ignore]
fn simil_0_1(structure: &str) -> Vec<BitArraySolution> {
    let lattice = Lattice::from_dft_json(
        format!("../exports/scaling_lattices/{}", structure),
        1.1,
        true,
    );

    let options = bit_array_settings!(
        lattice,
        solve_filter = BitArrayFilter::Flipped,
        difference_distance = 0.1,
        max_singlets = 0
    );

    let bit_lattice = lattice.get_intermediary(options);

    let result = bit_lattice.solve(true, true);
    print!(
        "\r{}{}",
        cursor::Right(19),
        clear::AfterCursor
    );
    result
}

#[divan::bench(max_time = 120, args = ["T04.json", "T08.json", "T09.json", "T12.json", "T12alt.json", "T16.json", "T20.json", "T25.json"])]
// #[ignore]
fn simil_0_5(structure: &str) -> Vec<BitArraySolution> {
    let lattice = Lattice::from_dft_json(
        format!("../exports/scaling_lattices/{}", structure),
        1.1,
        true,
    );

    let options = bit_array_settings!(
        lattice,
        solve_filter = BitArrayFilter::Flipped,
        difference_distance = 0.5,
        max_singlets = 0
    );

    let bit_lattice = lattice.get_intermediary(options);

    let result = bit_lattice.solve(true, true);
    print!(
        "\r{}{}",
        cursor::Right(19),
        clear::AfterCursor
    );
    result
}

fn main() {
    divan::main();
}
