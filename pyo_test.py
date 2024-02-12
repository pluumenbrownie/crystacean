from lattice_solver_python import Lattice, test_module
from findthosepoints import full_lattice_from_basis_vectors
import matplotlib.pyplot as plt

test_module()

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

# lattice = Lattice(bigger_boundary_points)

# basis_vectors = [[1.5, 0.0, 0.0], [0.75, 1.299038105676658, 0.0], [0.0, 0.0, 1.5]]


custom_boundary_points = full_lattice_from_basis_vectors(size=2)

lattice = Lattice(custom_boundary_points)
plt.plot( *lattice.points_to_plot(), "o")
plt.show()
bit_lattice = lattice.get_intermediary()
print(bit_lattice)
exit()
solutions = bit_lattice.solve(True)

to_show = {0, 8, 9}
for number, solution in enumerate(solutions):
    solved_lattice = lattice.to_solved_lattice(solution)
    # solved_lattice.export("exports", f"test_{number:>04}.json")

    if number in to_show:
        print(f"Solution {number}")
        plt.plot( *solved_lattice.points_to_plot(), "o")
        for num, point in enumerate(zip( *lattice.points_to_plot())):
            plt.annotate(str(num), (point[0], point[1]))

        plt.plot( *solved_lattice.midpoints_to_plot(), "x")
        plt.plot( *solved_lattice.tripoints_to_plot(), "s")
        plt.plot( *solved_lattice.singlets_to_plot(), "^", markersize=10)
        for num, point in enumerate(zip( *solved_lattice.oxygens_to_plot())):
            plt.annotate(str(num), (point[0], point[1]))
        plt.show()
# print(rust_those_points.print_in_rust((1, 2)))
