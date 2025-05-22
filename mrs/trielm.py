# -*- coding: utf-8 -*-
"""Triangle element.
"""
import numpy as np
from triellipt import amr
from triellipt import trimesh


class TriElm:
    """Triangle element.
    """

    def __init__(self, verts, nodes):

        self.verts = np.atleast_1d(verts)  # points
        self.nodes = np.atleast_1d(nodes)  # triangs

        self.data = {}
        self.mesh = None
        self.rank = 0

    @classmethod
    def from_data(cls, verts, nodes):
        return cls(verts, nodes).with_mesh()

    @property
    def anchor(self):
        return self.verts[0].real, self.verts[0].imag

    def with_mesh(self):

        self.mesh = trimesh.TriMesh.from_data(
            self.verts, np.array([[0, 1, 2]])
        )

        return self

    def makedata(self, key, func):

        self.data[key] = func(
            *self.mesh.points2d
        )

        return self

    def refine(self):

        unit = amr.getunit(self.mesh).with_data(self.data)
        unit = unit.refine()

        self.mesh = unit.mesh.twin()
        self.data = unit.data.copy()

        self.rank = self.rank + 1
        return self.aligned()

    def coarsen(self):

        unit = amr.getunit(self.mesh).with_data(self.data)

        unit.mesh = unit.mesh.shuffled(
            unit.find_subset(unit.mesh.size, self.anchor)
        )

        unit = unit.coarsen(
            unit.mesh.supertriu().reduce().trinums
        )

        self.mesh = unit.mesh.twin()
        self.data = unit.data.copy()

        self.rank = self.rank - 1
        return self.aligned()

    def aligned(self):

        new_mesh = self.mesh.alignnodes(self.anchor)
        dat_perm = new_mesh.meta['nodes-permuter']

        for key in self.data:
            self.data[key] = self.data[key][dat_perm]

        self.mesh = new_mesh.twin()
        return self

    def get_edges(self, shift=0):
        """Numbers of nodes across three edges as a packed triple.
        """

        nodnums = np.arange(
            3 * (2 ** self.rank)  # for aligned elements only!
        )

        corners = np.arange(
            0, nodnums.size, 2**self.rank
        )

        repnums = np.ones_like(nodnums)
        repnums[corners] = 2

        nodnums = np.roll(
            np.repeat(nodnums, repnums), -1
        )

        return nodnums + shift
