# -*- coding: utf-8 -*-
"""Cyclic path.
"""
import numpy as np
from triellipt.geom import loop


class CycPath:
    """Cyclic polygonal path.

    Properties
    ----------

     Name      | Description
    -----------|----------------------
    `points`   | Nodes positions.
    `numbers`  | Nodes numbers.

    """

    NODE = np.dtype(
        [('i', np.int32), ('z', np.complex128)]
    )

    def __init__(self, nodes=None):
        self.nodes = nodes

    @classmethod
    def from_path(cls, path):
        """Creates a cycle from a polygonal path.
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
    def points2d(self):
        return _unpack_complex(self.points)

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

    def dissect(self, angle):
        """Splits the cycle based on rotation angle.

        Parameters
        ----------
        angle : float
            Threshold angle for a node to become a corner.

        Returns
        -------
        PathMap
            Partition of a cycle.

        """

        corners, = np.where(
            self.angles() >= angle
        )

        return self.split_at_corners(corners)

    def split(self, bins):
        """Splits the cycle based on bins.

        Parameters
        ----------
        bins : flat-int-array-like
            Seed values of splitting bins.

        Returns
        -------
        PathMap
            Partition of a cycle.

        """
        return self.split_at_corners(np.add.accumulate(bins))

    def split_at_corners(self, corners):

        if corners.size == 0:
            return loop.PathMap.from_paths(self.points)

        paths = np.split(
            self.points, corners
        )

        paths = [
            path for path in paths if path.size != 0
        ]

        return loop.PathMap.from_paths(*paths)


def _triplet_angles(zone, ztwo, ztri):

    zdisp = ztwo - zone

    return np.angle(
        np.abs(zdisp) * (ztri - ztwo) / zdisp
    )


def _unpack_complex(argz):
    return np.vstack(
        [argz.real, argz.imag]
    )
