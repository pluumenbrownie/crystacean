use crystacean_rs::{
    bit_array_settings,
    test_points::{huge_points, lattice_points},
    BitArrayFilter, BitArraySettings, BitArraySolution, Lattice,
};
use termion::{cursor, clear};

fn main() {
    for margin in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0].into_iter().rev() {
        let results = similarity(margin);
        println!("\rMargin: {:<4} --- Structures: {}{}", margin, results.len(), clear::AfterCursor);
    }
    {
        let results = none();
        println!("\rMargin: {:<4} --- Structures: {}{}", "None", results.len(), clear::AfterCursor);
    }
}

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

    bit_lattice.solve(true, false)
}

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

    bit_lattice.solve(true, false)
}