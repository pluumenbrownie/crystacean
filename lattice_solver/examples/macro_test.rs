use lattice_solver::{BitArraySettings, Lattice, bit_array_settings, BitArrayFilter, SettingsBuilder};

fn main() {
    let lattice_points = vec![
        (
            vec![1.5, 0.0], 
            vec![vec![1.5, 5.196152422706632]]
        ),
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
        (
            vec![4.5, 0.0], 
            vec![vec![4.5, 5.196152422706632]]
        ),
        (
            vec![3.0, 0.0], 
            vec![vec![3.0, 5.196152422706632]]
        ),
        (vec![5.25, 1.299038105676658], vec![]),
        (vec![3.75, 1.299038105676658], vec![]),
        (vec![4.5, 2.598076211353316], vec![]),
        (vec![3.0, 2.598076211353316], vec![]),
        (vec![5.25, 3.897114317029974], vec![]),
        (vec![3.75, 3.897114317029974], vec![]),
    ];
    let lattice = Lattice::python_new(lattice_points, 1.1, true);
    let filter = bit_array_settings!(
        lattice,
        max_singlets = 1
    );

    println!("{filter:#?}");
}