use criterion::{criterion_group, criterion_main, Criterion};
use lattice_solver::{bit_array_settings, test_points::lattice_points, BitArraySettings, Lattice};
use std::time::Duration;

fn criterion_long(c: &mut Criterion) {
    let lattice = Lattice::python_new(lattice_points(), 1.1, true);
    let options = bit_array_settings!(&lattice);
    let bit_lattice = lattice.get_intermediary(options);
    let mut test_bitarray = bit_lattice.get_bitarray();
    test_bitarray.set(0, true);
    // test_bitarray.set(0, true);

    c.bench_function("full_run", |b| b.iter(|| bit_lattice.solve(true, false)));
    // c.bench_function("fib 20", |b| b.iter(|| fibonacci(black_box(20))));
}

criterion_group!(
    name=full_run;
    config=Criterion::default()
    .warm_up_time(Duration::from_secs(1))
    .measurement_time(Duration::from_secs(180));
    targets=criterion_long
);

criterion_main!(full_run);
