use crystacean_rs::{bit_array_settings, BitArrayFilter, BitArraySolution, Lattice};
use termion::{clear, cursor};



#[divan::bench(max_time = 120, args = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])]
fn similarity(margin: f32) -> Vec<BitArraySolution> {
    let lattice = Lattice::from_dft_json(
        "../test_lattices/T20.json".into(),
        1.1,
        true,
    );

    let options = bit_array_settings!(
        lattice,
        solve_filter = BitArrayFilter::Flipped,
        difference_distance = margin,
        max_singlets = 0
    );

    let bit_lattice = lattice.get_intermediary(options);

    let result = bit_lattice.solve(true, true);
    print!(
        "\r{}{}",
        cursor::Right(15),
        clear::AfterCursor
    );
    result
}


#[divan::bench(max_time = 120)]
fn none() -> Vec<BitArraySolution> {
    let lattice = Lattice::from_dft_json(
        "../test_lattices/T20.json".into(),
        1.1,
        true,
    );

    let options = bit_array_settings!(
        lattice,
        solve_filter = BitArrayFilter::None,
        max_singlets = 0
    );

    let bit_lattice = lattice.get_intermediary(options);

    let result = bit_lattice.solve(true, true);
    print!(
        "\r{}{}",
        cursor::Right(15),
        clear::AfterCursor
    );
    result
}

fn main() {
    divan::main();
}
