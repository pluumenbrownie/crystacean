from crystacean import Lattice, test_module
from basis_vectors import full_lattice_from_basis_vectors
import matplotlib.pyplot as plt
from matplotlib import rcParams
from tqdm import tqdm

rcParams.update({'font.size': 11})
lattice_size = 10
site_size = 16

test_module()

custom_boundary_points = full_lattice_from_basis_vectors(size=1)
lattice = Lattice(custom_boundary_points, 1.1)

# for num, point in enumerate(zip( *lattice.points_to_plot())):
#     plt.annotate(str(num), (point[0], point[1]))

plt.plot( *lattice.points_to_plot(), "o", markersize=lattice_size, label="Previous layer")
plt.plot( *lattice.midpoints_to_plot(), "x", markersize=site_size, label="Midpoint", markeredgewidth=4)
plt.plot( *lattice.tripoints_to_plot(), "s", markersize=site_size, label="Tripoint")
plt.plot( *lattice.singlets_to_plot(), "^", markersize=site_size, label="Single")
plt.axis('equal')
# plt.xlim(-0.5, 8.5)
# plt.ylim(-0.5, 6)
plt.xlabel("x (Å)")
plt.ylabel("y (Å)")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.11),
          ncol=4, fancybox=True)
plt.show()

bit_lattice = lattice.get_intermediary()
print(bit_lattice)

solutions = bit_lattice.solve(True)

progress = tqdm(
    enumerate(solutions),
    total=len(solutions),
    bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}  "
)

for number, solution in progress:
    solved_lattice = lattice.to_solved_lattice(solution)
    # solved_lattice.export("exports/T6_noloops", f"T6_{number:>06}.json")
    # solved_lattice.export_as_ase_json(f"exports/T6_noloops/T6_{number:>06}.json")

    # if number in to_show:
    # solved_lattice.export_as_ase_json(f"exports/sizeable_remade/sizeable_rust_{number:>06}.json")
    progress.set_description(desc=f"Solution {number}")
#     # print(f"Solution {number}")
    plt.plot( *solved_lattice.points_to_plot(), "o", markersize=lattice_size, label="Previous layer")
#     for num, point in enumerate(zip( *lattice.points_to_plot())):
#         plt.annotate(str(num), (point[0], point[1]))

    plt.plot( *solved_lattice.midpoints_to_plot(), "x", markersize=site_size, label="Midpoint", markeredgewidth=4)
    plt.plot( *solved_lattice.tripoints_to_plot(), "s", markersize=site_size, label="Tripoint")
    plt.plot( *solved_lattice.singlets_to_plot(), "^",  markersize=site_size, label="Single")
#     for num, point in enumerate(zip( *solved_lattice.oxygens_to_plot())):
#         plt.annotate(str(num), (point[0], point[1]))
    plt.axis('equal')
    # plt.xlim(-0.5, 8.5)
    # plt.ylim(-0.5, 6)
    plt.xlabel("x (Å)")
    plt.ylabel("y (Å)")
    # plt.title("")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.11),
          ncol=4, fancybox=True)
    plt.show()
