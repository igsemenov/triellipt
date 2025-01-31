# -*- coding: utf-8 -*-
"""Nodal masses.
"""
import numpy as np
from triellipt.amr import utils_


def getmasser(mesh):
    return TriMasser.from_mesh(mesh)


class TriMasser:
    """Evaluator of nodal masses.

    Attributes
    ----------
    mesh : TriMesh
        Parent triangle mesh.
    meta : dict
        Masser metadata.

    """

    def __init__(self, mesh, meta):
        self.mesh = mesh
        self.meta = meta

    @classmethod
    def from_mesh(cls, mesh):

        nodesmap = mesh.nodesmap()
        triareas = utils_.mesh_areas(mesh)

        meta = {
            'nodesmap': nodesmap,
            'triareas': triareas[nodesmap.trinums]
        }

        return cls(mesh, meta)

    @property
    def count(self):
        return self.mesh.npoints

    @property
    def nnodes(self):
        return self.nodesmap.data.nitems

    @property
    def triareas(self):
        return self.meta['triareas']

    @property
    def nodesmap(self):
        return self.meta['nodesmap']

    @property
    def nodes_range(self):
        return self.nodesmap.nodes_range.astype(int)

    def areas(self):
        """Computes areas of the nodes.
        """
        return self.mass_diag(np.ones(self.count))

    def mass_full(self, data):
        """Computes masses using the full mass-matrix.
        """
        _ = data[self.mesh.triangs]
        return self.mass_full_from_data_table(_)

    def mass_full_from_data_table(self, data_table):
        mass12, mass13 = self.get_mass_on_edges(data_table)
        return self.sum_over_nodesmap(mass12 + mass13)

    def get_mass_on_edges(self, data_table):

        data_self = data_table.flat[self.stream_self]
        data_ccw1 = data_table.flat[self.stream_ccw1]
        data_ccw2 = data_table.flat[self.stream_ccw2]

        data12 = 0.5 * (data_self + data_ccw1)
        data13 = 0.5 * (data_self + data_ccw2)

        mass12 = data12 * (self.triareas / 6.)
        mass13 = data13 * (self.triareas / 6.)

        return mass12, mass13

    def mass_diag(self, data):
        """Computes masses using the lumped mass-matrix.
        """
        _ = data[self.mesh.triangs]
        return self.mass_diag_from_data_table(_)

    def mass_diag_from_data_table(self, data_table):
        mass_nodes = self.get_mass_at_nodes(data_table)
        return self.sum_over_nodesmap(mass_nodes)

    def get_mass_at_nodes(self, data_table):

        data_nodes = data_table.flat[self.stream_self]
        mass_nodes = data_nodes * (self.triareas / 3.)

        return mass_nodes

    def sum_over_nodesmap(self, data_long):
        return np.add.reduceat(
            data_long, self.bins_reduce
        )

    @property
    def bins_reduce(self):
        return self.nodesmap.data.bins_reduce

    @property
    def stream_self(self):
        return 3 * self.nodesmap.trinums + self.nodesmap.locnums

    @property
    def stream_ccw1(self):
        return 3 * self.nodesmap.trinums + self.nodesmap.locnums1

    @property
    def stream_ccw2(self):
        return 3 * self.nodesmap.trinums + self.nodesmap.locnums2
