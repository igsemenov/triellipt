# -*- coding: utf-8 -*-
"""FEM unit partition.
"""
import itertools as itr
import numpy as np
from triellipt.fem import femvector, femmatrix


def getpartt(unit, spec):
    """Creates a FEM unit partition.

    Parameters
    ----------
    spec : dict
        Partition specification.

    Returns
    -------
    FEMPartt
        Partition object.

    """
    partter = UnitPartter.from_unit(unit)
    return partter.getpartt(spec)


def make_base_partition(unit):
    """Makes the edge-core unit partition.
    """

    loops = [
        loop.nodnums_unique for loop in unit.loops
    ]

    meta = {
        'name': 'base',
        'dirichlet-sides': [
            i + 1 for i in range(len(loops))
        ]
    }

    edges = {
        i + 1: loop for (i, loop) in enumerate(loops)
    }

    return FEMPartt(
        unit, meta, edges
    )


class UnitAgent:
    """Operator on a FEM unit.
    """

    def __init__(self, unit):
        self.unit = unit
        self.cache = {}

    @classmethod
    def from_unit(cls, unit):
        return cls(unit)

    def get_loops_nodenums(self):

        nodenums = [
            loop.nodnums1 for loop in self.unit.loops
        ]

        return np.hstack(nodenums)

    def get_loops_points(self):

        points = [
            loop.nodes_complex for loop in self.unit.loops
        ]

        return np.hstack(points)


class UnitPartter(UnitAgent):
    """Makes unit partition.
    """

    def getpartt(self, spec):

        self.cache['spec'] = spec

        return FEMPartt(
            self.unit, spec, self.make_edges()
        )

    def make_edges(self):

        sides = self.split_loops()

        return {
            i + 1: v for i, v in enumerate(sides)
        }

    def split_loops(self):

        sides = itr.chain.from_iterable(
            self.split_loop(l) for l in self.loops_as_arrays
        )

        return sides

    def split_loop(self, loop):

        anchors = self.find_anchors()

        _, inds, _ = np.intersect1d(
            loop, anchors, assume_unique=True, return_indices=True
        )

        return np.split(loop, inds)

    def find_anchors(self):

        points = self.get_loops_points()
        nodenums = self.get_loops_nodenums()

        closestinds = [
            _find_closest(points, p) for p in self.anchors_complex
        ]

        return nodenums[closestinds]

    @property
    def loops_as_arrays(self):
        return [
            loop.nodnums1 for loop in self.unit.loops
        ]

    @property
    def anchors_complex(self):
        return [
            complex(x, y) for x, y in self.cache['spec']['anchors']
        ]

    @property
    def anchors_nodnums(self):
        return self.cache['anchors-nodnums']


FEMParttError = type(
    'FEMParttError', (Exception,), {}
)


class DataPartt:
    """Partition data.

    Attributes
    ----------
    unit : FEMUnit
        Parent FEM unit.
    meta : dict
        Partition metadata.
    edge : dict
        Edge sides numbered from one.

    """

    def __init__(self, unit, meta, edge):
        self.unit = unit
        self.meta = meta
        self.edge = edge

    def __getitem__(self, key):
        if key == 0:
            return self.core
        if key in self.edge:
            return self.edge[key]
        raise FEMParttError(
            f"no edge {key} in the partition '{self.name}'"
        )

    @property
    def name(self):
        return self.meta['name']

    @property
    def dirich_sides(self):
        return self.meta['dirichlet-sides']

    @property
    def core(self):
        return self.getcore()

    def getcore(self):

        if not self.dirich_sides:
            return np.arange(
                self.unit.mesh_count
            )

        nodes_to_remove = [
            self.edge[ind] for ind in self.dirich_sides
        ]

        return np.delete(
            np.arange(self.unit.mesh_count), np.hstack(nodes_to_remove)
        )

    def nodes2d(self, key):
        return self.unit.mesh.points2d[:, self[key]]

    def nodes_complex(self, key):
        return self.unit.mesh.points[self[key]]


class FEMPartt(DataPartt):
    """FEM unit partition.

    Properties
    ----------

    Name        | Description
    ------------|-------------------------
    `core`      | Core section.
    `edge`      | Map of edge sections.
    `meta`      | Partition metadata.

    """

    def new_vector(self):
        """Creates a new FEM vector.

        Returns
        -------
        VectorFEM
            New empty FEM vector.

        """
        return femvector.getvector(self)

    def new_matrix(self, operator, add_constr):
        """Creates a new FEM matrix.

        Parameters
        ----------
        operator : flat-float-array
            Linear combination of the basic FEM operators.
        add_constr : bool
            Constraints are included in the matrix, if True.

        Returns
        -------
        MatrixFEM
            Resulting FEM matrix.

        """
        if add_constr is False:
            return self.make_free_matrix(operator)
        return self.make_full_matrix(operator)

    def make_free_matrix(self, data):
        body, meta = self.unit.factory_free.feed_data(data)
        return femmatrix.getmatrix(self, body, meta)

    def make_full_matrix(self, data):
        body, meta = self.unit.factory_full.feed_data(data)
        return femmatrix.getmatrix(self, body, meta)

    def get_nodes(self, key):
        """Retrieves the points of the partition section.

        Parameters
        ----------
        key : int
            Number of the partition section.

        Returns
        -------
        two-row-float-array
            Points of the partition section stacked horizontally.

        """
        return self.unit.mesh.points2d[:, self[key]]


def _find_closest(points, anchor):
    return np.argmin(abs(points - anchor))
