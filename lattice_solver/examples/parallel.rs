use crystacean_rs::{bit_array_settings, test_points::lattice_points, BitArraySettings, Lattice};

fn main() {
    let lattice = Lattice::python_new(lattice_points(), 1.1, true);
    let options = bit_array_settings!(lattice);
    let bit_lattice = lattice.get_intermediary(options);

    let _ = bit_lattice.solve_parallel(true, false);
}
