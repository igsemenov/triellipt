# -*- coding: utf-8 -*-
"""Cyclic path.
"""
import numpy as np


class CycPath:
    """Cyclic path with focus on geometry.

    Properties
    ----------

     Name      | Description
    -----------|---------------------------------------------
    `points`   | Nodes positions.
    `numbers`  | Nodes numbers.
    `previnds` | Cyclic index with the shift −1
    `nextinds` | Cyclic index with the shift +1

    """

    NODE = np.dtype(
        [('i', np.int32), ('z', np.complex128)]
    )

    def __init__(self, nodes=None):
        self.nodes = nodes

    @classmethod
    def from_path(cls, path):
        """Creates a cycle from a polygonal chain.
        """

        nodes = np.array(
            list(enumerate(path)), dtype=cls.NODE
        )

        return cls(nodes)

    @property
    def size(self):
        return len(self.nodes)

    @property
    def numbers(self):
        """Global nodes numbers.
        """
        return self.nodes[:]['i']

    @property
    def points(self):
        """Path points.
        """
        return self.nodes[:]['z']

    @property
    def prev_points(self):
        """Left side neighbours of path points.
        """
        return self.nodes[self.previnds]['z']

    @property
    def next_points(self):
        """Right side neighbours of path points.
        """
        return self.nodes[self.nextinds]['z']

    @property
    def previnds(self):
        return np.r_[-1, np.arange(self.size-1)]

    @property
    def nextinds(self):
        return np.r_[np.arange(1, self.size), 0]

    def angles(self):
        """Provides rotation angles of edges.

        Returns
        -------
        flat-float-array
            Rotation angles of edges at each node.

        """
        return _triplet_angles(
            self.prev_points, self.points, self.next_points
        )


def _triplet_angles(zone, ztwo, ztri):

    zdisp = ztwo - zone

    return np.angle(
        np.abs(zdisp) * (ztri - ztwo) / zdisp
    )
