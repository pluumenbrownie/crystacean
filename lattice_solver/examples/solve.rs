#![allow(clippy::excessive_precision)]

use lattice_solver::{
    bit_array_settings, BitArrayFilter, BitArraySettings, Lattice
};

fn main() {
    let lattice_points = vec![
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
    ];
    // let lattice_points = vec![
    //     (
    //         vec![3.076, 0.0], vec![vec![3.076, 5.327788284081867]]
    //     ),
    //     (
    //         vec![0.0, 0.0], vec![
    //             vec![0.0, 5.327788284081867],
    //             vec![6.152, 0.0],
    //             vec![6.152, 5.327788284081867]
    //         ]
    //     ),
    //     (vec![4.614, 2.6638941420409337], vec![]),
    //     (
    //         vec![1.538, 2.6638941420409337], vec![vec![7.69, 2.6638941420409337]]
    //     )
    // ];
    let lattice = Lattice::python_new(lattice_points, 1.1, true);

    let options = bit_array_settings!(
        lattice
    );
    let bit_lattice = lattice.get_intermediary(options);

    // bit_lattice.print_distances();
    println!("{}", bit_lattice.__str__());

    let solutions = bit_lattice.solve(true, false);
    println!("Solutions found: {}", solutions.len());

    // let solutions_filtered = bit_lattice.solve_filtered(true, false);
    // println!("Solutions found: {}", solutions.len());
    // for (number, solution) in solutions.iter().enumerate() {
    //     let solution_lattice = lattice.to_solved_lattice(solution);
    //     solution_lattice.export(&OsString::from("../exports/test4"), format!("test2_{number}.json"));
    // }
}
