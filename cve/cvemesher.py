# -*- coding: utf-8 -*-
"""CVE mesh maker.
"""
import numpy as np
from triellipt.utils import tables


def getcvemesh(mesh):
    maker = CVEMeshMaker.from_mesh(mesh)
    return maker.make_cvemesh()


class CVEMeshMaker:
    """Maker of the CVE mesh.
    """

    def __init__(self, mesh=None):
        self.mesh = mesh
        self.cache = {}

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    def make_cvemesh(self):

        _ = self.make_cvemesh_data()
        _ = self.push_cvemesh()

        return _

    def make_cvemesh_data(self):

        _ = self.make_cve_triangs()
        _ = self.make_cve_points()

        return True

    def push_cvemesh(self):

        cvemesh = self.mesh.from_data(
            self.cache['cve-points'],  self.cache['cve-triangs']
        )

        cvemesh.meta['edgescodes'] = self.cache['edgespaired'].copy()
        return cvemesh

    def make_cve_triangs(self):

        edges = self.mesh.edges_paired()

        self.cache['edgespaired'] = edges
        self.cache['cve-triangs'] = tables.norm_table(edges)

        return True

    def make_cve_points(self):

        _ = self.make_edgesnodes()
        _ = self.flat_edgesnodes()

        return _

    def make_edgesnodes(self):

        edges = 0.5 * (
            self.trinodes1 + self.trinodes2
        )

        self.cache['edges-nodes'] = edges
        return True

    def flat_edgesnodes(self):

        points = np.zeros(
            np.amax(self.cache['cve-triangs']) + 1, dtype=complex
        )

        nodes = self.cache['edges-nodes']
        triangs = self.cache['cve-triangs']

        points[triangs.flat] = nodes.flat
        self.cache['cve-points'] = points.copy()

        return True

    @property
    def trinodes1(self):
        return self.mesh.points[
            self.mesh.triangs[:, [0, 1, 2]]
        ]

    @property
    def trinodes2(self):
        return self.mesh.points[
            self.mesh.triangs[:, [1, 2, 0]]
        ]
