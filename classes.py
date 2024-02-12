from __future__ import annotations
from itertools import chain
import json
from typing import Sequence
from dataclasses import dataclass, field
import matplotlib.pylab as plt
import numpy as np


SINGLE_POINT_ENERGY = 1.4
MID_POINT_ENERGY = 0.7
TRI_POINT_ENERGY = 0.4


@dataclass
class LatticePoint:
    x: float
    y: float
    _connected_to: None | SinglePoint | MidPoint | TriPoint = field(default=None, init=False)

    def get_location(self) -> list[float]:
        return [self.x, self.y]
    
    def get_real_location(self) -> list[float]:
        return self.get_location()

    def get_link(self) -> LatticePoint:
        return self
    
    @property
    def connected_to(self) -> None | SinglePoint | MidPoint | TriPoint:
        return self._connected_to
    
    def set_connection(self, connection: None | SinglePoint | MidPoint | TriPoint):
        if not self._connected_to and not connection:
            raise ValueError("Connection already None")
        if self._connected_to and connection:
            raise ValueError("Point already connected")
        self._connected_to = connection

    
    def __hash__(self) -> int:
        return hash((self.x, self.y, self.connected_to))
    

@dataclass
class GhostPoint(LatticePoint):
    linked_point: LatticePoint

    def get_link(self) -> LatticePoint:
        return self.linked_point
    
    @property
    def connected_to(self) -> None | SinglePoint | MidPoint | TriPoint:
        return self.linked_point.connected_to

    def set_connection(self, connection: None | SinglePoint | MidPoint | TriPoint):
        self.linked_point.set_connection(connection)
    
    def __hash__(self) -> int:
        return hash((self.x, self.y, self.connected_to))


@dataclass
class SinglePoint:
    x: float
    y: float
    _connections: list[LatticePoint] = field(init=False)
    populated: bool = False

    def set_connections(self, connections: list[LatticePoint]) -> None:
        assert len(connections) == 1
        self._connections = connections

    def get_connections(self) -> list[LatticePoint]:
        return self._connections
    
    def available(self) -> bool:
        if self.populated:
            return False
        
        for connection in self.get_connections():
            if connection.connected_to:
                return False
        return True
    
    def populate(self) -> bool:
        """
        Connect oxygen to connections, if possible. 
        Returns `True` if all connections are available, else returns `False` and 
        doesn't connect.
        """
        assert not self.populated, "This point is already populated."

        for connection in self.get_connections():
            if connection.connected_to:
                return False
        
        for connection in self.get_connections():
            connection.set_connection(self)
        self.populated = True
        return True
    
    def depopulate(self) -> None:
        """
        Reset the point and connected lattice sites.
        """
        assert self.populated, "Point not populated."

        for connection in self.get_connections():
            connection.set_connection(None)
        self.populated = False

    def __hash__(self) -> int:
        return hash((self.x, self.y))


@dataclass
class MidPoint(SinglePoint):
    def set_connections(self, connections: list[LatticePoint]) -> None:
        assert len(connections) == 2
        self._connections = connections

    def __hash__(self) -> int:
        return hash((self.x, self.y))


@dataclass
class TriPoint(SinglePoint):
    def set_connections(self, connections: list[LatticePoint]) -> None:
        assert len(connections) == 3
        self._connections = connections

    def __hash__(self) -> int:
        return hash((self.x, self.y))


def points_to_plot(input_points: Sequence[Plottable]) -> tuple[list[float], list[float]]:
    """
    `Plottable = LatticePoint | SinglePoint | MidPoint | TriPoint`
    """
    x_points = [point.x for point in input_points]
    y_points = [point.y for point in input_points]
    return x_points, y_points


@dataclass
class FullLattice:
    points: list[LatticePoint]
    singles: list[SinglePoint]
    midpoints: list[MidPoint]
    tripoints: list[TriPoint]


    def energy(self) -> float:
        total = 0
        for tri in self.tripoints:
            if tri.populated:
                total += TRI_POINT_ENERGY
        for mid in self.midpoints:
            if mid.populated:
                total += MID_POINT_ENERGY
        for single in self.singles:
            if single.populated:
                total += SINGLE_POINT_ENERGY
        return total
    
    def is_solved(self) -> bool:
        for point in self.points:
            if not point.connected_to:
                return False
        return True
    
    def empty_amount(self) -> int:
        total = 0
        for point in self.points:
            if not point.connected_to:
                total += 1
        return total
    
    def plot(self, draw_single: bool = False):
        plt.plot( *points_to_plot(self.points), "o")
        for num, point in enumerate(self.points):
            plt.annotate(str(num), (point.x, point.y))

        plt.plot( *points_to_plot(self.midpoints), "x")
        plt.plot( *points_to_plot(self.tripoints), "s")
        if draw_single:
            plt.plot( *points_to_plot(self.singles), "^", markersize=10) 
        plt.show()
    
    def plot_ghost_connections(self):
        plt.plot( *points_to_plot(self.points), "o")
        for num, point in enumerate(self.points):
            plt.annotate(str(num), (point.x, point.y))

        for point in self.points:
            linked_point = point.get_link()
            if not linked_point == point:
                plt.plot( *points_to_plot([linked_point, point]), "--")
        plt.show()
    
    def as_spg_tuple(self, oxygens: bool = True) -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float]], list[int]]:
        basis_vectors = [(0.0, 0.0, 1.5)]
        return_points = []
        return_elements = []

        for silicon in self.points:
            if isinstance(silicon, GhostPoint):
                print("skipped")
                continue
            return_points.append(( *silicon.get_location(), 0.0))
            return_elements.append(14)
            if len(basis_vectors) == 1:
                if not silicon.y == 0.0:
                    basis_vectors.append((silicon.x, silicon.y, 0.0))
            elif len(basis_vectors) == 2:
                if silicon.y == 0.0:
                    basis_vectors.append((silicon.x, silicon.y, 0.0))

        if oxygens:
            for oxygen in chain(self.tripoints, self.midpoints, self.singles):
                return_points.append((oxygen.x, oxygen.y, 1.5))
                return_elements.append(8)
        
        return basis_vectors, return_points, return_elements

    def export(self, filename: str):
        export_dict = {
            "lattice_points" : [],
            "tripoints": [],
            "midpoints": [],
            "singles": []
        }
        for point in self.points:
            export_dict["lattice_points"].append(
                {
                    "x": point.x,
                    "y": point.y,
                    "ghost": isinstance(point, GhostPoint)
                }
            )
        for tri in self.tripoints:
            export_dict["tripoints"].append(
                {
                    "x": tri.x,
                    "y": tri.y
                }
            )
        for mid in self.midpoints:
            export_dict["midpoints"].append(
                {
                    "x": mid.x,
                    "y": mid.y
                }            )
        for single in self.singles:
            export_dict["singles"].append(
                {
                    "x": single.x,
                    "y": single.y
                }
            )
        with open(f"exports/{filename}.json", mode="w") as file:
            json.dump(export_dict, file, indent=4)
        
    def __hash__(self) -> int:
        return hash(tuple(self.points))


def from_file(filename: str) -> FullLattice:
    points: list[LatticePoint] = []
    singles: list[SinglePoint] = []
    midpoints: list[MidPoint] = []
    tripoints: list[TriPoint] = []
    with open(filename) as file:
        structure = json.load(file)
    for point in structure["lattice_points"]:
        points.append(LatticePoint(point["x"], point["y"]))
    for single in structure["singles"]:
        singles.append(SinglePoint(single["x"], single["y"], populated=True))
    for mid in structure["midpoints"]:
        midpoints.append(MidPoint(mid["x"], mid["y"], populated=True))
    for tri in structure["tripoints"]:
        tripoints.append(TriPoint(tri["x"], tri["y"], populated=True))
    return FullLattice(points, singles, midpoints, tripoints)


type Plottable = LatticePoint | SinglePoint | MidPoint | TriPoint

if __name__ == '__main__':
    import os
    path = "exports/sizeable_unique_1/(4, 0, 4)"
    for structure in os.listdir(path):
        lattice = from_file(f"{path}/{structure}")
        lattice.plot(draw_single=True)
    # lattice = from_file("exports/sizeable_unique_wrong/(5, 0, 1)/sizeable_rust_0022.json")
    # lattice.plot(draw_single=True)