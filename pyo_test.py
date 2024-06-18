from crystacean import Lattice, test_module, from_dft_json
from basis_vectors import full_lattice_from_basis_vectors
import matplotlib.pyplot as plt
from matplotlib import rcParams
from tqdm import tqdm

rcParams.update({'font.size': 16})
lattice_size = 10
site_size = 16

test_module()

tiny_points = [
    ([0.0, 0.0], []),
    ([1.6, 2.7], []),
    ([3.1, 0.0], []),
    ([4.6, 2.7], []),
]

boundary_points = [
    ([0.0, 0.0], [
        [0.0, 5.3],
        [6.2, 5.3],
        [6.2, 0.0],
    ]),
    ([1.6, 2.7], [
        [7.7, 2.7],
    ]),
    ([3.1, 0.0], [
        [3.1, 5.3],
    ]),
    ([4.6, 2.7], []),
]


bigger_boundary_points = [
    ([0.0, 0.0], [
        [0.0, 10.7],
        [9.3, 0.0],
        [9.3, 10.7],
    ]),
    ([0.0, 5.3], [
        [9.3, 5.3],
    ]),
    ([1.6, 2.7], [
        [10.8, 2.7],
    ]),
    ([1.6, 8.0], [
        [10.8, 8.0],
    ]),
    ([3.1, 0.0], [
        [3.1, 10.7],
    ]),
    ([3.1, 5.3], []),
    ([4.6, 2.7], []),
    ([4.6, 8.0], []),
    ([6.2, 0.0], [
        [6.2, 10.7],
    ]),
    ([6.2, 5.3], []),
    ([7.7, 2.7], []),
    ([7.7, 8.0], []),
]

# lattice = Lattice(tiny_points, 1.1)

# basis_vectors = [[1.5, 0.0, 0.0], [0.75, 1.299038105676658, 0.0], [0.0, 0.0, 1.5]]


# custom_boundary_points = full_lattice_from_basis_vectors(size=2)
# lattice = Lattice(custom_boundary_points, 1.1)

lattice = from_dft_json("DFT_results/T1S.json", 3.5, False)
# lattice.diagnostic_ase()
# for num, point in enumerate(zip( *lattice.points_to_plot())):
#     plt.annotate(str(num), (point[0], point[1]))

# plt.plot( *lattice.points_to_plot(), "o", markersize=lattice_size)
# plt.plot( *lattice.midpoints_to_plot(), "x", markersize=site_size)
# plt.plot( *lattice.tripoints_to_plot(), "s", markersize=site_size)
# plt.plot( *lattice.singlets_to_plot(), "^", markersize=site_size)
# plt.axis('equal')
# plt.xlim(0, 5)
# plt.ylim(-5, 3)
# plt.xlabel("x (Å)")
# plt.ylabel("y (Å)")
# # for num, x, y in lattice.no_rings_plot():
# #     plt.annotate(str(num), (x, y))
# plt.show()

# noloops = lattice.no_rings()

bit_lattice = lattice.get_intermediary()
print(bit_lattice)
# filtered_bit = bit_lattice.filtered(noloops)
# print(filtered_bit)
# exit()
solutions = bit_lattice.solve(True)
# solutions = filtered_bit.solve(True)
exit()
progress = tqdm(
    enumerate(solutions),
    total=len(solutions),
    bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}  "
)

# to_show = {0, 19}
for number, solution in progress:
    solved_lattice = lattice.to_solved_lattice(solution)
    # solved_lattice.export("exports/T6_noloops", f"T6_{number:>06}.json")
    # solved_lattice.export_as_ase_json(f"exports/T6_noloops/T6_{number:>06}.json")

    # if number in to_show:
    # solved_lattice.export_as_ase_json(f"exports/sizeable_remade/sizeable_rust_{number:>06}.json")
    progress.set_description(desc=f"Solution {number}")
#     # print(f"Solution {number}")
    plt.plot( *solved_lattice.points_to_plot(), "o", markersize=lattice_size)
#     for num, point in enumerate(zip( *lattice.points_to_plot())):
#         plt.annotate(str(num), (point[0], point[1]))

    plt.plot( *solved_lattice.midpoints_to_plot(), "x", markersize=site_size)
    plt.plot( *solved_lattice.tripoints_to_plot(), "s", markersize=site_size)
    plt.plot( *solved_lattice.singlets_to_plot(), "^",  markersize=site_size)
#     for num, point in enumerate(zip( *solved_lattice.oxygens_to_plot())):
#         plt.annotate(str(num), (point[0], point[1]))
    plt.axis('equal')
    plt.xlim(0, 5)
    plt.ylim(-5, 3)
    plt.xlabel("x (Å)")
    plt.ylabel("y (Å)")
    # plt.title("")
    # plt.legend(loc=)
    # plt.show()
    
    plt.show()
# print(rust_those_points.print_in_rust((1, 2)))
