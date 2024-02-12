use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;
use lattice_solver::Lattice;

fn criterion_benchmark(c: &mut Criterion) {
    let lattice_points = vec![(vec![1.5, 0.0], vec![vec![1.5, 5.196152422706632]]), (vec![0.0, 0.0], vec![vec![0.0, 5.196152422706632], vec![6.0, 0.0], vec![6.0, 5.196152422706632]]), (vec![2.25, 1.299038105676658], vec![]), (vec![0.75, 1.299038105676658], vec![vec![6.75, 1.299038105676658]]), (vec![1.5, 2.598076211353316], vec![]), (vec![0.0, 2.598076211353316], vec![vec![6.0, 2.598076211353316]]), (vec![2.25, 3.897114317029974], vec![]), (vec![0.75, 3.897114317029974], vec![vec![6.75, 3.897114317029974]]), (vec![4.5, 0.0], vec![vec![4.5, 5.196152422706632]]), (vec![3.0, 0.0], vec![vec![3.0, 5.196152422706632]]), (vec![5.25, 1.299038105676658], vec![]), (vec![3.75, 1.299038105676658], vec![]), (vec![4.5, 2.598076211353316], vec![]), (vec![3.0, 2.598076211353316], vec![]), (vec![5.25, 3.897114317029974], vec![]), (vec![3.75, 3.897114317029974], vec![])];
    let lattice = Lattice::python_new(lattice_points);
    let bit_lattice = lattice.get_intermediary();
    let mut test_bitarray = bit_lattice.get_bitarray();
    test_bitarray.set(0, true);
    // test_bitarray.set(0, true);
    
    c.bench_function("get_possibilities", |b| b.iter(|| bit_lattice.get_possibilities(&test_bitarray)));
    // c.bench_function("fib 20", |b| b.iter(|| fibonacci(black_box(20))));
    
}

criterion_group!(
name=benches; 
config=Criterion::default()
    .warm_up_time(Duration::from_secs(1))
    .measurement_time(Duration::from_secs(1));
targets=criterion_benchmark
);

criterion_main!(benches);