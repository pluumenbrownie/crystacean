use crystacean_rs::{BitArraySettings, Lattice, test_points::lattice_points, bit_array_settings, BitArrayFilter};

fn main() {
    let lattice = Lattice::python_new(lattice_points(), 1.1, true);
    let filter = bit_array_settings!(
        lattice,
        max_singlets = 1
    );
    let correct = BitArraySettings::create(1, 0.05, lattice.find_max(), BitArrayFilter::None);

    assert_eq!(filter, correct);
}