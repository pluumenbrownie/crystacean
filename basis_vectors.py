import numpy as np

from itertools import product


def full_lattice_from_basis_vectors(size: int = 1) -> list[tuple[list[float], list[list[float]]]]:
    points = []
    bottom_points = {}
    left_points = {}
    vec_x = np.array([1.5, 0.0, 0.0]) / 1.5 * 3.076
    vec_y = np.array([0.75, 1.299038105676658, 0.0]) / 1.5 * 3.076
    vec_z = np.array([0.0, 0.0, 1.5])

    for x, y in product(range(0, size), range(0, size)):
        # create main cluster
        new_point_lb = (list((2*x - y)     * vec_x + 2*y      * vec_y)[0:2], [])
        new_point_rb = (list((2*x - y  + 1)* vec_x + 2*y      * vec_y)[0:2], [])
        new_point_lt = (list((2*x - y)     * vec_x +(2*y + 1) * vec_y)[0:2], [])
        new_point_rt = (list((2*x - y  + 1)* vec_x +(2*y + 1) * vec_y)[0:2], [])
        points.append(new_point_rb)
        points.append(new_point_lb)
        points.append(new_point_rt)
        points.append(new_point_lt)

        # note left/bottom boundaries
        if y == 0:
            bottom_points[round(new_point_lb[0][0], 4)] = new_point_lb
            bottom_points[round(new_point_rb[0][0], 4)] = new_point_rb
        if x == 0:
            left_points[round(new_point_lb[0][1], 4)] = new_point_lb
            left_points[round(new_point_lt[0][1], 4)] = new_point_lt

        # create and link top/right boundaries
        if y == size - 1:
            bound_point_1 = list((2*x - y - 1) * vec_x + 2*(y + 1)*vec_y)[0:2]
            bount_point_2 = list((2*x - y) * vec_x + 2*(y + 1)* vec_y)[0:2]
            bottom_points[round(bound_point_1[0], 4)][1].append(bound_point_1)
            bottom_points[round(bount_point_2[0], 4)][1].append(bount_point_2)
        if x == size - 1:
            bound_point_3 = list((2*x - y  + 2)* vec_x + 2*y      * vec_y)[0:2]
            bount_point_4 = list((2*x - y  + 2)* vec_x +(2*y + 1) * vec_y)[0:2]
            left_points[round(bound_point_3[1], 4)][1].append(bound_point_3)
            left_points[round(bount_point_4[1], 4)][1].append(bount_point_4)
            if y == size - 1:
                corner_point = list((2*x - y  + 1)* vec_x +(2*y + 2) * vec_y)[0:2]
                bottom_points[0.0][1].append(corner_point)

    return points