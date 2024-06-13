#![allow(clippy::excessive_precision)]
use lattice_solver::{BitArraySettings, Lattice};
use termion::{clear, cursor};

fn lattice_points() -> Vec<(Vec<f32>, Vec<Vec<f32>>)> {
    vec![
        (vec![1.5, 0.0], vec![vec![1.5, 5.196152422706632]]),
        (
            vec![0.0, 0.0],
            vec![
                vec![0.0, 5.196152422706632],
                vec![6.0, 0.0],
                vec![6.0, 5.196152422706632],
            ],
        ),
        (vec![2.25, 1.299038105676658], vec![]),
        (
            vec![0.75, 1.299038105676658],
            vec![vec![6.75, 1.299038105676658]],
        ),
        (vec![1.5, 2.598076211353316], vec![]),
        (
            vec![0.0, 2.598076211353316],
            vec![vec![6.0, 2.598076211353316]],
        ),
        (vec![2.25, 3.897114317029974], vec![]),
        (
            vec![0.75, 3.897114317029974],
            vec![vec![6.75, 3.897114317029974]],
        ),
        (vec![4.5, 0.0], vec![vec![4.5, 5.196152422706632]]),
        (vec![3.0, 0.0], vec![vec![3.0, 5.196152422706632]]),
        (vec![5.25, 1.299038105676658], vec![]),
        (vec![3.75, 1.299038105676658], vec![]),
        (vec![4.5, 2.598076211353316], vec![]),
        (vec![3.0, 2.598076211353316], vec![]),
        (vec![5.25, 3.897114317029974], vec![]),
        (vec![3.75, 3.897114317029974], vec![]),
    ]
}

#[divan::bench(max_time = 120, args = [true, false])]
fn full_run(silent: bool) -> Vec<lattice_solver::BitArraySolution> {
    let lattice_points = lattice_points();

    let lattice = Lattice::python_new(lattice_points, 1.1, true);
    let options = BitArraySettings::create(2, 0.1, lattice.find_max());

    let bit_lattice = lattice.get_intermediary(options);

    let result = bit_lattice.solve(true, silent);
    print!("\r{}{}", cursor::Right(22), clear::AfterCursor);
    result
}

#[divan::bench(max_time = 120, args = [true, false])]
fn full_run_parallel(silent: bool) -> Vec<lattice_solver::BitArraySolution> {
    let lattice_points = lattice_points();

    let lattice = Lattice::python_new(lattice_points, 1.1, true);

    let options = BitArraySettings::create(2, 0.1, lattice.find_max());
    let bit_lattice = lattice.get_intermediary(options);

    let result = bit_lattice.solve_parallel(true, silent);
    print!("\r{}{}", cursor::Right(22), clear::AfterCursor);
    result
}

fn main() {
    divan::main();
}
