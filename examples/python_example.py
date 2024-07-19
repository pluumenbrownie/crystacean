import os
from crystacean import from_dft_json

# Make the output directory
save_to = "exports/T16_example"
os.mkdir(save_to)
# Import the lattice from the file
lattice = from_dft_json("test_lattices/T16.json", 3.5, False)
# get_intermediary creates a solver class and can be supplied with a few options
bit_lattice = lattice.get_intermediary(
    # no singlets
    max_singlets = 0,
    # use the similatiry filter, with a margin of 0.1 angstrom
    use_filter = True,
    difference_distance = 0.1,
)
# Create a list bitvectors representing valid solutions
# The solve method gets passed True to find all solutions, instead of stopping
# after finding the first one
solutions = bit_lattice.solve(True)
# Save all of the structures
for number, solution in enumerate(solutions):
    # Solutions need to be reverted back to lattices for exporting
    solved_lattice = lattice.to_solved_lattice(solution)
    # Save the files as json
    solved_lattice.export_as_ase_json(
        f"example_{number:>04}.json", save_to
    )