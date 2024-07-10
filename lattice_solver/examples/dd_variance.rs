use lattice_solver::{
    bit_array_settings,
    test_points::{huge_points, lattice_points},
    BitArrayFilter, BitArraySettings, BitArraySolution, Lattice,
};
use termion::{cursor, clear};

fn main() {
    for margin in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9].into_iter().rev() {
        let results = similarity(margin);
        println!("\rMargin: {:<4} --- Structures: {}{}", margin, results.len(), clear::AfterCursor);
    }
}

fn similarity(margin: f32) -> Vec<BitArraySolution> {
    let lattice = Lattice::from_dft_json(
        "../exports/scaling_lattices/T20.json".into(),
        1.1,
        true,
    );

    let options = bit_array_settings!(
        lattice,
        solve_filter = BitArrayFilter::Similarity,
        difference_distance = margin,
        max_singlets = 0
    );

    let bit_lattice = lattice.get_intermediary(options);

    bit_lattice.solve(true, false)
}