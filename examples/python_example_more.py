import os
import matplotlib.pyplot as plt
from crystacean import from_dft_json

os.mkdir("exports/T16_example")
lattice = from_dft_json("test_lattices/T16.json", 3.5, False)
bit_lattice = lattice.get_intermediary(
    max_singlets = 0,
    use_filter = True,
    difference_distance = 0.1,
)

# To apply the no rings filter:
noloops = lattice.no_rings()
bit_lattice = bit_lattice.filtered(noloops)

solutions = bit_lattice.solve(True)
for number, solution in enumerate(solutions):
    solved_lattice = lattice.to_solved_lattice(solution)
    solved_lattice.export_as_ase_json(
        f"example_{number:>04}.json", save_to
    )

    # Structures can be plotted with matplotlib as follows
    plt.plot( *solved_lattice.points_to_plot())
    plt.plot( *solved_lattice.midpoints_to_plot(), "x")
    plt.plot( *solved_lattice.tripoints_to_plot(), "s")
    plt.plot( *solved_lattice.singlets_to_plot(), "^")
    plt.axis("equal")
    plt.show()