use crystacean_rs::{bit_array_settings, BitArrayFilter, Lattice};

fn main() {
    // Read a lattice from a file with good standard values.
    let lattice = Lattice::from_dft_json("../test_lattices/T16.json".into(), 1.1, true);

    // The options for the bit array are given in a macro, to allow for default values
    let options = bit_array_settings!(
        // First argument should be the `Lattice`
        lattice,
        max_singlets = 0,
        // BitArrayFilter::Flipped is the fastest filter, and the one 
        // used in the python module
        solve_filter = BitArrayFilter::Flipped,
        difference_distance = 0.1
    );

    let bit_lattice = lattice.get_intermediary(options);
    // The first argument in solve determines if the function will stop after finding
    // the first solution, the second argument can hide the progress bar.
    let solutions = bit_lattice.solve(true, false);
    
    for (number, solution) in solutions.iter().enumerate() {
        let solved_lattice = lattice.to_solved_lattice(solution);
        solved_lattice.export_as_ase_json(
            (&format!("../exports/T16_example/example_{number:>4}.json"))
        );
    }
}