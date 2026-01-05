# -*- coding: utf-8 -*-
"""FEM-DtN unit.
"""
import numpy as np
from triellipt.fem import femunit


def getdtn(mesh, anchors, mode=None):
    """Creates a FEM-DtN unit.

    Parameters
    ----------
    mesh : TriMesh
        Input triangle mesh.
    anchors : Iterable
        Provides `(float, float)` points to partition the mesh boundary.
    mode : str = None
        Solver mode — "fvm" or "fem" (default).

    Returns
    -------
    FEMDtN
        FEM-DtN computing unit.

    Notes
    -----

    FEM-DtN units are only supported for simply-connected meshes.

    """

    if len(anchors) == 0:
        raise FEMDtNError(
            'anchors is empty. Provide at least one anchor.'
        )

    if mode is None:
        return getdtn_fem(mesh, anchors)

    if mode == "fvm":
        return getdtn_fvm(mesh, anchors)
    return getdtn_fem(mesh, anchors)


def getdtn_fem(mesh, anchors):

    maker = FEMDtNMaker.from_mesh(mesh)
    maker = maker.with_mode('fem').with_anchors(anchors)

    return maker.getdtn()


def getdtn_fvm(mesh, anchors):

    maker = FEMDtNMaker.from_mesh(mesh)
    maker = maker.with_mode('fvm').with_anchors(anchors)

    return maker.getdtn()


class FEMDtN:
    """FEM-DtN computing unit.

    Properties
    ----------

    Name        | Description
    ------------|----------------------
    `fem`       | Parent FEM unit.
    `dtn`       | DtN partition.

    """

    def __init__(self, unit=None):
        self.unit = unit

    @property
    def fem(self):
        return self.unit

    @property
    def dtn(self):
        return self.unit.partts['dtn']

    @property
    def dirich_sides_(self):
        return self.dirich_sides()

    def dirich_sides(self) -> list:
        """Returns a list of Dirichlet sides.
        """
        return self.dtn.meta['dirichlet-sides']

    def switch_side(self, key):
        """Switches the DtN state of the boundary part.

        Parameters
        ----------
        key : int
            Index of the boundary part.

        Returns
        -------
        self
            Unit itself.

        """

        if key == 0:
            return self

        if key not in self.dirich_sides_:
            self.dirich_sides_.append(key)
            return self

        self.dirich_sides_.remove(key)
        return self


class MeshAgent:
    """Operator on a mesh.
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self.mode = None
        self.meta = {}

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    def with_mode(self, mode):
        self.mode = mode
        return self

    def with_anchors(self, anchors):
        self.meta['anchors'] = anchors
        return self

    @property
    def anchors(self):
        return self.meta['anchors']

    @property
    def nanchors(self):
        return len(self.anchors)


class FEMDtNMaker(MeshAgent):
    """Makes a FEM-DtN unit from a mesh. 
    """

    def getdtn(self):

        unit = self.make_fem_unit()
        unit = self.add_partition(unit)

        return FEMDtN(unit)

    def make_fem_unit(self):

        unitmaker = self.get_unitmaker()

        unit = unitmaker(
            self.mesh, anchors=[self.anchors[0]]
        )

        if len(unit.loops) != 1:
            raise FEMDtNError(
                'cannot make a FEMDtN unit, mesh is not simply-connected'
            )

        loop = unit.loops[0]

        self.meta['loop2d'] = loop.nodes2d
        self.meta['loop_z'] = loop.nodes_complex

        return unit

    def get_unitmaker(self):
        if self.mode == "fvm":
            return femunit.getunit_fvm
        return femunit.getunit_fem

    def add_partition(self, unit):

        unit = unit.add_partition(
            self.get_partition_spec()
        )

        return unit

    def get_partition_spec(self):

        anchors = self.get_partition_anchors()

        drsides = list(
            np.arange(2 * self.nanchors) + 1
        )

        return {
            'name': 'dtn',
            'anchors': anchors,
            'dirichlet-sides': drsides
        }

    def get_partition_anchors(self):

        indices = self.get_partition_indices()
        corners = self.get_partition_corners(indices)

        anchors = [
            corners[0]
        ]

        for ind in range(1, self.nanchors):
            anchors += [
                self.anchors[ind], corners[ind]
            ]

        return anchors

    def get_partition_indices(self):

        dists = [
            np.abs(self.loop_z - complex(*p)) for p in self.anchors
        ]

        return [
            np.argmin(v) for v in dists
        ]

    def get_partition_corners(self, indices):
        return [
            tuple(self.loop2d[:, i + 1]) for i in indices
        ]

    @property
    def loop2d(self):
        return self.meta['loop2d']

    @property
    def loop_z(self):
        return self.meta['loop_z']


FEMDtNError = type(
    'FEMDtNError', (Exception,), {}
)
