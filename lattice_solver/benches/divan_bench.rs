#![allow(clippy::excessive_precision)]
use lattice_solver::{
    Lattice, BitArraySolution
};
use termion::{clear, cursor};

pub fn lattice_points() -> Vec<(Vec<f32>, Vec<Vec<f32>>)> {
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

#[divan::bench(max_time = 120)]
fn huge_run() -> Vec<BitArraySolution> {
    println!();
    let lattice_points = lattice_points();
    
    let lattice = Lattice::python_new(lattice_points, 1.1, true);

    let bit_lattice = lattice.get_intermediary();

    let result = bit_lattice.solve(true);
    print!("\r{}{}{}", cursor::Up(1), cursor::Right(18), clear::AfterCursor);
    result
}

fn main() {
    divan::main();
}
