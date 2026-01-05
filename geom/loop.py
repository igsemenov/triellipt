# -*- coding: utf-8 -*-
"""Loop of curves.
"""
import itertools as itr
import numpy as np
from triellipt.geom import dump

CurvesLoopError = type(
    'CurvesLoopError', (Exception,), {}
)


class CurvesLoop:
    """Loop of connected curves.

    Attributes
    ----------
    curves : tuple
        Sequence of connected curves.

    Properties
    ----------

     Name         | Description
    --------------|----------------------
    `startpoints` | Curves start points.
    `endpoints`   | Curves end points.

    """

    TOL = 1e-12

    def __init__(self, curves=None):
        self.curves = curves

    @property
    def ncurves(self):
        return len(self.curves)

    @property
    def startpoints(self):
        return [
            curve.startpoint for curve in self.curves
        ]

    @property
    def endpoints(self):
        return [
            curve.endpoint for curve in self.curves
        ]

    @classmethod
    def from_curves(cls, *curves):

        if len(curves) == 0:
            raise CurvesLoopError(
                "expected at least one curve, got 0"
            )

        loop = cls(curves)

        if not loop.is_closed():
            raise CurvesLoopError(
                "loop is not closed, check the curves"
            )

        index = loop.find_break()

        if index is None:
            return loop

        curv0 = index
        curv1 = index + 1

        raise CurvesLoopError(
            f"curves {curv0} and {curv1} are disconnected"
        )

    def is_closed(self):
        return self.are_curves_linked(
            self.curves[-1], self.curves[0]
        )

    def find_break(self):

        contacts = itr.pairwise(self.curves)

        for index, curves in enumerate(contacts):
            if not self.are_curves_linked(*curves):
                return index
        return None

    @classmethod
    def are_points_equal(cls, point1, point2):
        return np.isclose(
            point1, point2, atol=cls.TOL, rtol=0.0
        )

    @classmethod
    def are_curves_linked(cls, curve1, curve2):
        return cls.are_points_equal(
            curve1.endpoint, curve2.startpoint
        )

    def discretize(self, *params):
        """Discretizes the curves loop.

        Parameters
        ----------
        params : *(nparts, ratio)
            Partition parameters for each curve, see `Curve().partition()`

        Returns
        -------
        PathsMap
            Polygonal path colored at curves.

        """

        if len(params) == 0:
            return np.array([], dtype=complex)

        args = self.zip_partition_args(*params)

        paths = [
            curve.partition(*params) for curve, params in args
        ]

        paths = [
            path[:-1] for path in paths
        ]

        return PathMap.from_paths(*paths)

    def zip_partition_args(self, *params):

        nparams = len(params)

        if nparams == 0:
            return None
        if nparams >= self.ncurves:
            return zip(self.curves, params)
        return itr.zip_longest(
            self.curves, params, fillvalue=params[-1]
        )


class PathMap:
    """Polygonal path with colored nodes.

    Properties
    ----------

     Name      | Description
    -----------|------------------------------
    `colors`   | Nodes colors.
    `numbers`  | Nodes numbers.
    `points`   | Nodes positions (complex).
    `points2d` | Nodes positions (xy-rows).

    """

    NODE = np.dtype(
        [('i', np.int32), ('c', np.int32), ('z', np.complex128)]
    )

    def __init__(self, nodes=None):
        self.nodes = nodes

    @property
    def size(self):
        return len(self.nodes)

    @classmethod
    def from_paths(cls, *paths):

        if not paths:
            return cls()

        counts = [
            len(path) for path in paths
        ]

        nums = np.arange(sum(counts))

        keys = np.repeat(
            np.arange(len(paths)), counts
        )

        data = np.hstack(paths)

        nodes = np.array(
            list(zip(nums, keys, data)), dtype=cls.NODE
        )

        return cls(nodes)

    @property
    def numbers(self):
        return self.nodes[:]['i']

    @property
    def colors(self):
        return self.nodes[:]['c']

    @property
    def points(self):
        return self.nodes[:]['z']

    @property
    def points2d(self):
        return np.vstack(
            [self.points.real, self.points.imag]
        )

    def getcopy(self):
        return self.__class__(np.copy(self.nodes))

    def togeo(self, geopath, seeds) -> None:
        """Dumps path to the geo file.

        Parameters
        ----------
        geopath : str
            Absolute path to the `.geo` file.
        seeds : dict
            Maps colours to the seed mesh sizes.

        """
        return dump.dump_path_map(self, geopath, seeds)

    def atcolors(self, *colors):
        """Fetches a subpath with the specified colors.

        Parameters
        ----------
        colors : *int
            Colors in the subpath.

        Returns
        -------
        PathMap
            The resulting subpath.

        """

        mask = np.isin(
            self.colors, colors
        )

        new_nodes = self.nodes[mask]
        return self.__class__(new_nodes)

    def repaint(self, color, newcolor):
        """Changes the specified color to the new one.
        """
        self.colors[self.colors == color] = newcolor
        return self.getcopy()

    def rshift(self, color1, color2):
        """Shifts the contact of two colors to the right by one node.
        """

        contact = self.find_contact(color1, color2)

        if contact is None:
            return self.getcopy()

        loc1, loc2 = contact

        self.colors[loc2] = self.colors[loc1]
        return self.getcopy()

    def lshift(self, color1, color2):
        """Shifts the contact of two colors to the left by one node.
        """

        contact = self.find_contact(color1, color2)

        if contact is None:
            return self.getcopy()

        loc1, loc2 = contact

        self.colors[loc1] = self.colors[loc2]
        return self.getcopy()

    def find_contact(self, color1, color2):
        return self.cyccolors.find_contact(color1, color2)

    @property
    def cyccolors(self):
        return _CycleColor.from_colors(self.colors)


class _CycleColor:

    COLOR = np.dtype(
        [('i', np.int32), ('c', np.int32)]
    )

    def __init__(self, data=None):
        self.data = data

    @classmethod
    def from_colors(cls, colors):

        data = np.array(
            list(enumerate(colors)), dtype=cls.COLOR
        )

        cycdata = np.r_[data, data[0]]
        return cls(cycdata)

    @property
    def cycinds(self):
        return self.data[:]['i']

    @property
    def colors(self):
        return self.data[:]['c']

    def find_contact(self, color1, color2):

        locs1, locs2 = self.get_contacts()

        match1 = self.colors[locs1] == color1
        match2 = self.colors[locs2] == color2

        touchs = np.logical_and(match1, match2)

        if not any(touchs):
            return None

        touch = np.argmax(touchs)

        loc1 = locs1[touch]
        loc2 = locs2[touch]

        return self.cycinds[loc1], self.cycinds[loc2]

    def get_contacts(self):

        locis = list(
            self.gen_contacts()
        )

        return np.array(locis).T

    def gen_contacts(self):
        for dat1, dat2 in itr.pairwise(self.data):
            if dat1['c'] != dat2['c']:
                yield dat1['i'], dat2['i']
