from __future__ import annotations
from copy import deepcopy
from typing import Sequence
import matplotlib.pyplot as plt
from matplotlib import rcParams
from itertools import combinations, permutations
from scipy.spatial import KDTree
from dataclasses import dataclass, field
from classes import *

from basis_vectors import full_lattice_from_basis_vectors


def find_ox_sites(lattice: list[list[float]]) -> FullLattice:
    lattice.sort(key=lambda x: 100*x[0] + x[1])

    lattice_points = [LatticePoint(point[0], point[1]) for point in lattice]
    single_points = [SinglePoint(point.x, point.y) for point in lattice_points]
    for s_point, point in zip(single_points, lattice_points):
        s_point.set_connections([point])

    point_kdtree = KDTree(lattice) # type: ignore
    node_search_distance = point_kdtree.query(lattice_points[0].get_location(), k=[2])[0][0] * 1.1
    unsorted_pairs: set[tuple[int, int]] = point_kdtree.query_pairs(node_search_distance, output_type='set')
    pairs = sorted(unsorted_pairs)
    print(node_search_distance)

    midpoints: list[MidPoint] = []
    for pair in pairs:
        point1 = lattice[pair[0]]
        point2 = lattice[pair[1]]
        new_point = MidPoint((point1[0] + point2[0]) / 2, (point1[1] + point2[1]) / 2)
        new_point.set_connections([
            lattice_points[pair[0]], 
            lattice_points[pair[1]]
        ])
        midpoints.append(new_point)

    tripoints: list[TriPoint] = []
    for pair1, pair2 in combinations(pairs, 2):
        if (not pair1[0] == pair2[0]):
            continue
        point1 = lattice[pair1[0]]
        point2 = lattice[pair1[1]]
        point3 = lattice[pair2[1]]
        # if abs(point2[0] - point3[0]) > node_search_distance or abs(point2[1] - point3[1]) > node_search_distance:
        if not (pair1[1], pair2[1]) in unsorted_pairs:
            continue

        # print(f"{pair1=}, {pair2};      {point1=}, {point2=}, {point3=}")
        new_point = TriPoint((point1[0] + point2[0] + point3[0]) / 3, (point1[1] + point2[1]+ point3[1]) / 3)
        new_point.set_connections( [
            lattice_points[pair1[0]], 
            lattice_points[pair1[1]], 
            lattice_points[pair2[1]]
        ] )
        tripoints.append(new_point)

    return FullLattice(lattice_points, single_points, midpoints, tripoints)


def point_sorting_key(point: LatticePoint) -> float:
    return point.x * 100 + point.y


def find_ox_sites_bounds(lattice: list[tuple[list[float], list[list[float]]]]) -> FullLattice:
    """
    Where `ghost_lattice` is a superset of `lattice`.
    """
    lattice.sort(key=lambda x: 100*x[0][0] + x[0][1])

    lattice_points: list[LatticePoint] = []
    for point in lattice:
        point_object = LatticePoint(point[0][0], point[0][1])
        lattice_points.append(point_object)
        for ghost in point[1]:
            lattice_points.append(GhostPoint(ghost[0], ghost[1], point_object))


    single_points = [SinglePoint(point.x, point.y) for point in lattice_points if not isinstance(point, GhostPoint)]
    for s_point, point in zip(single_points, [point for point in lattice_points if not isinstance(point, GhostPoint)]):
        s_point.set_connections([point.get_link()])
    
    lattice_points.sort(key=point_sorting_key)
    new_lattice: list[list[float]] = []
    for point in lattice:
        new_lattice.append(point[0])
        for ghost in point[1]:
            new_lattice.append(ghost)
    new_lattice.sort(key=lambda x: 100*x[0] + x[1])

    for one, other in zip(lattice_points, new_lattice):
        assert one.x == other[0]
        assert one.y == other[1]

    point_kdtree = KDTree(new_lattice) # type: ignore
    node_search_distance = point_kdtree.query(lattice_points[0].get_location(), k=[2])[0][0] * 1.1

    unsorted_pairs: set[tuple[int, int]] = point_kdtree.query_pairs(node_search_distance, output_type='set')
    pairs = sorted(unsorted_pairs)
    print(node_search_distance)

    midpoints: list[MidPoint] = []
    for pair in pairs:
        point1 = lattice_points[pair[0]]
        point2 = lattice_points[pair[1]]
        if isinstance(point1, GhostPoint) and isinstance(point2, GhostPoint):
            continue
        new_point = MidPoint((point1.x + point2.x) / 2, (point1.y + point2.y) / 2)
        new_point.set_connections([
            point1.get_link(),
            point2.get_link()
        ])
        midpoints.append(new_point)

    tripoints: list[TriPoint] = []
    for pair1, pair2 in combinations(pairs, 2):
        if (not pair1[0] == pair2[0]):
            continue
        if not (pair1[1], pair2[1]) in unsorted_pairs:
            continue

        point1 = lattice_points[pair1[0]]
        point2 = lattice_points[pair1[1]]
        point3 = lattice_points[pair2[1]]
        # if abs(point2[0] - point3[0]) > node_search_distance or abs(point2[1] - point3[1]) > node_search_distance:
        if isinstance(point1, GhostPoint) and isinstance(point2, GhostPoint) and isinstance(point2, GhostPoint):
            continue

        # print(f"{pair1=}, {pair2};      {point1=}, {point2=}, {point3=}")
        new_point = TriPoint((point1.x + point2.x + point3.x) / 3, (point1.y + point2.y+ point3.y) / 3)
        new_point.set_connections( [
            point1.get_link(), 
            point2.get_link(), 
            point3.get_link()
        ] )
        tripoints.append(new_point)

    return FullLattice(lattice_points, single_points, midpoints, tripoints)


def points_to_plot(input_points: Sequence[Plottable]) -> tuple[list[float], list[float]]:
    """
    `Plottable = LatticePoint | SinglePoint | MidPoint | TriPoint`
    """
    x_points = [point.x for point in input_points]
    y_points = [point.y for point in input_points]
    return x_points, y_points


def make_site_plot(
    lattice: FullLattice,
    draw_single: bool = False
):
    plt.plot( *points_to_plot(lattice.points), "o")
    for num, point in enumerate(lattice.points):
        plt.annotate(str(num), (point.x, point.y))

    plt.plot( *points_to_plot(lattice.midpoints), "x")
    plt.plot( *points_to_plot(lattice.tripoints), "s")
    if draw_single:
        plt.plot( *points_to_plot(lattice.singles), "^", markersize=10) 
    plt.show()


def fill_oxygen_maxtri(lattice: FullLattice) -> FullLattice:
    """
    Fill the lattice with points by first trying to populate all `TriPoint`s.
    """
    # shuffle(tripoints)
    filled_tri: list[TriPoint] = []
    for tri in lattice.tripoints:
        if tri.populate():
            filled_tri.append(tri)
    
    filled_mid: list[MidPoint] = []
    for mid in lattice.midpoints:
        if mid.populate():
            filled_mid.append(mid)

    filled_single: list[SinglePoint] = []
    for single in lattice.singles:
        if single.populate():
            filled_single.append(single)
    
    return FullLattice(lattice.points, filled_single, filled_mid, filled_tri)


def fill_oxygen_midonly(lattice: FullLattice) -> FullLattice:
    """
    Fill the lattice with points by first trying to populate all `MidPoint`s.
    """
    # shuffle(midpoints)
    filled_mid: list[MidPoint] = []
    for mid in lattice.midpoints:
        if mid.populate():
            filled_mid.append(mid)

    filled_single: list[SinglePoint] = []
    for single in lattice.singles:
        if single.populate():
            filled_single.append(single)
    
    return FullLattice(lattice.points, filled_single, filled_mid, [])


@dataclass
class Solver:
    lattice: FullLattice
    lowest_energy = 1_000_000.0
    archive: set[int] = field(default_factory=set, init=False)
    depth: int = field(default=0, init=False)

    def start_solve(self, depth: int = 0) -> list[FullLattice]:
        current_gen: list[FullLattice] = []
        current_gen.append(deepcopy(self.lattice))
        gen_nr = 0
        solutions: list[FullLattice] = []
        while not solutions or (depth and gen_nr <= depth):
            if depth:
                print(f"Generation {gen_nr}/{depth}")
            else:
                print(f"Generation {gen_nr}")

            next_gen = []
            for lattice in current_gen:
                options_found = 0
                spare_lattice = deepcopy(lattice)
                for number, tripoint in enumerate(lattice.tripoints):
                    if not tripoint.available():
                        continue
                    spare_lattice.tripoints[number].populate()
                    spare_lattice_hash = hash(spare_lattice)
                    if spare_lattice_hash in self.archive or spare_lattice.energy() >= self.lowest_energy:
                        spare_lattice.tripoints[number].depopulate()
                        continue

                    self.archive.add(spare_lattice_hash)
                    if spare_lattice.is_solved():
                        solutions.append(deepcopy(spare_lattice))
                    elif not solutions or (depth and gen_nr < depth):
                        next_gen.append(deepcopy(spare_lattice))
                    options_found += 0
                    spare_lattice.tripoints[number].depopulate()

                for number, midpoint in enumerate(lattice.midpoints):
                    if not midpoint.available():
                        continue
                    spare_lattice.midpoints[number].populate()
                    spare_lattice_hash = hash(spare_lattice)
                    if spare_lattice_hash in self.archive or spare_lattice.energy() >= self.lowest_energy:
                        spare_lattice.midpoints[number].depopulate()
                        continue

                    self.archive.add(spare_lattice_hash)
                    if spare_lattice.is_solved():
                        solutions.append(deepcopy(spare_lattice))
                    elif not solutions or (depth and gen_nr < depth):
                        next_gen.append(deepcopy(spare_lattice))
                    options_found += 0
                    spare_lattice.midpoints[number].depopulate()

                # do not try to fit singles in if not necessary
                if options_found:
                    continue
                for number, single in enumerate(lattice.singles):
                    if not single.available():
                        continue
                    spare_lattice.singles[number].populate()
                    spare_lattice_hash = hash(spare_lattice)
                    if spare_lattice_hash in self.archive or spare_lattice.energy() >= self.lowest_energy:
                        spare_lattice.singles[number].depopulate()
                        continue

                    self.archive.add(spare_lattice_hash)
                    if spare_lattice.is_solved():
                        solutions.append(deepcopy(spare_lattice))
                    elif not solutions or (depth and gen_nr < depth):
                        next_gen.append(deepcopy(spare_lattice))
                    spare_lattice.singles[number].depopulate()
            current_gen = next_gen
            gen_nr += 1
        
        if not solutions:
            raise ValueError("Depth was to shallow. Deepen search or remove depth requirement.")
        return solutions

    # def start_solve_accelerated(self, depth: int = 0) -> list[FullLattice]:
        


def reduce_lattice(lattice: FullLattice) -> FullLattice:
    return FullLattice(
        lattice.points,
        [single for single in lattice.singles if single.populated],
        [mid for mid in lattice.midpoints if mid.populated],
        [tri for tri in lattice.tripoints if tri.populated],
    )


if __name__ == '__main__':
    litterally_just_a_line = [
        [0.0, 0.0],
        [1.5, 0.0]
    ]
    equilateral_triangle_points = [
        [0.0, 0.0],
        [0.75, 1.299038105676658],
        [1.5, 0.0]
    ]
    flipped_triangle = [
        [0.0, 1.299038105676658],
        [0.75, 0.0],
        [1.5, 1.299038105676658]
    ]
    double_triangle_points = [
        [0.0, 0.0],
        [0.75, 1.299038105676658],
        [1.5, 0.0],
        [2.25, 1.299038105676658]
    ]
    big_lattice_points = [
        [0.0, 5.3],
        [1.6, 13.4],
        [3.1, 5.3],
        [9.3, 0.0],
        [10.8, 2.7],
        [7.7, 2.7],
        [1.6, 8.0],
        [6.2, 0.0],
        [13.9, 8.0],
        [12.4, 5.3],
        [3.1, 10.7],
        [0.0, 10.7],
        [4.6, 8.0],
        [13.9, 2.7],
        [4.6, 2.7],
        [3.1, 0.0],
        [10.8, 8.0],
        [9.3, 5.3],
        [10.8, 13.4],
        [6.2, 5.3],
        [1.6, 2.7],
        [0.0, 0.0],
        [6.2, 10.7],
        [13.9, 13.4],
        [12.4, 10.7],
        [7.7, 13.4],
        [7.7, 8.0],
        [12.4, 0.0],
        [4.6, 13.4],
        [9.3, 10.7]
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
    test_points = [
        (
            [1.5, 0.0], 
            [[1.5, 5.196152422706632]]
        ),
        (
            [0.0, 0.0],
            [
                [0.0, 5.196152422706632],
                [6.0, 0.0],
                [6.0, 5.196152422706632],
            ],
        ),
        ([2.25, 1.299038105676658], []),
        (
            [0.75, 1.299038105676658],
            [[6.75, 1.299038105676658]],
        ),
        ([1.5, 2.598076211353316], []),
        (
            [0.0, 2.598076211353316],
            [[6.0, 2.598076211353316]],
        ),
        ([2.25, 3.897114317029974], []),
        (
            [0.75, 3.897114317029974],
            [[6.75, 3.897114317029974]],
        ),
        (
            [4.5, 0.0], 
            [[4.5, 5.196152422706632]]
        ),
        (
            [3.0, 0.0], 
            [[3.0, 5.196152422706632]]
        ),
        ([5.25, 1.299038105676658], []),
        ([3.75, 1.299038105676658], []),
        ([4.5, 2.598076211353316], []),
        ([3.0, 2.598076211353316], []),
        ([5.25, 3.897114317029974], []),
        ([3.75, 3.897114317029974], []),
    ]
    test = full_lattice_from_basis_vectors(3)
    big_lattice_points = [point for point in big_lattice_points if point[0] < 12 and point[1] < 9]
    rcParams.update({'font.size': 11})
    lattice_size = 10
    site_size = 16

    # filled_lattice = find_ox_sites(big_lattice_points)
    filled_lattice = find_ox_sites_bounds(test)
    filled_lattice.plot_report(lattice_size, site_size)
    # points = full_lattice_from_basis_vectors(size=2)
    # filled_lattice = find_ox_sites_bounds(points)
    filled_lattice.plot_ghost_connections()
    # make_site_plot(filled_lattice, draw_single=True)

    # solved_lattice = fill_oxygen_maxtri(filled_lattice)
    # solved_lattice = fill_oxygen_midonly(filled_lattice)
    # shuffle(filled_lattice.midpoints)
    # shuffle(filled_lattice.tripoints)
    exit(0)
    solution_list = Solver(filled_lattice).start_solve(depth=0)
    
    # Use the "just look at it" theorum to select for unique solutions
    # exportset = {0, 1, 2, 3, 25, 26, 37, 70}
    # exportset = {0, 9, 10}
    for number, lattice in enumerate(solution_list):
        print(f"Showing {number+1} of {len(solution_list)}")
        print(f"Energy: {lattice.energy()}")

        reduced = reduce_lattice(lattice)
        reduced.plot_report(lattice_size, site_size)
        # if number in exportset:
        #     reduced.export(f"bigger_lattice_{number:0>4}")
