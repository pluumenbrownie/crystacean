from basis_vectors import full_lattice_from_basis_vectors
import spglib
from classes import *
from findthosepoints import find_ox_sites_bounds, find_ox_sites


def show_symmetry(symmetry):
    for i in range(symmetry["rotations"].shape[0]):
        print("  --------------- %4d ---------------" % (i + 1))
        rot = symmetry["rotations"][i]
        trans = symmetry["translations"][i]
        print("  rotation:")
        for x in rot:
            print("     [%2d %2d %2d]" % (x[0], x[1], x[2]))
        print("  translation:")
        print("     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2]))
        print(f"     ({trans[0]:8.5f} {trans[1]:8.5f} {trans[2]:8.5f})" )


def show_lattice(lattice):
    print("Basis vectors:")
    for vec, axis in zip(lattice, ("a", "b", "c")):
        print(
            "%s %10.5f %10.5f %10.5f"
            % (
                tuple(
                    axis,
                )
                + tuple(vec)
            )
        )


rutile = (
    [(4, 0, 0), (0, 4, 0), (0, 0, 3)],
    [
        (0, 0, 0),
        (0.5, 0.5, 0.5),
        (0.3, 0.3, 0.0),
        (0.7, 0.7, 0.0),
        (0.2, 0.8, 0.5),
        (0.8, 0.2, 0.5),
    ],
    [14, 14, 8, 8, 8, 8],
)

equilateral_triangle_sym = (
    [[0.0, 0.0, 1.5], [0.75, 1.299038105676658, 0.0], [1.5, 0.0, 0.0]],
    [[0.0, 0.0, 0.0], [0.75, 1.299038105676658, 0.0], [1.5, 0.0, 0.0]],
    [14, 14, 14]
)

driehoek = [
    [0.0, 0.0], 
    [0.75, 1.299038105676658], 
    [1.5, 0.0],
    # [2.25, 1.299038105676658], 
]

nonbound = [
    [0.0, 0.0],
    [1.6, 2.7],
    [3.1, 0.0],
    [4.6, 2.7]
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


# full_lattice = find_ox_sites(nonbound)
# full_lattice = find_ox_sites_bounds(boundary_points)
full_lattice = full_lattice_from_basis_vectors(size=2)
# spg = full_lattice.as_spg_tuple(oxygens=False)
spg = full_lattice.as_spg_tuple(oxygens=True)
symmetry = spglib.get_symmetry(spg)
full_lattice.plot(draw_single=True)
full_lattice.plot_ghost_connections()
print(symmetry)
show_symmetry(symmetry)
# print(full_lattice.as_spg_tuple())
# print("  Spacegroup of Rutile is %s." % spglib.get_spacegroup(rutile))
# print("  Pointgroup of Rutile is %s." % spglib.get_pointgroup(symmetry["rotations"])[0])
# show_lattice(rutile)