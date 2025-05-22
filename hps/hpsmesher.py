# -*- coding: utf-8 -*-
"""HPS mesh maker.
"""
import numpy as np
from triellipt.utils import tables


def gethpsmesh(mesh):
    maker = HPSMeshMaker.from_mesh(mesh)
    return maker.make_hpsmesh()


class HPSMeshMaker:
    """Maker of the HPS mesh.
    """

    def __init__(self, mesh=None):
        self.mesh = mesh
        self.cache = {}

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    def make_hpsmesh(self):

        _ = self.make_hpsmesh_data()
        _ = self.push_hpsmesh()

        return _

    def make_hpsmesh_data(self):

        triangs = self.make_hps_triangs()
        points = self.make_hps_points()

        return (triangs, points)

    def push_hpsmesh(self):

        hpsmesh = self.mesh.from_data(
            self.cache['hps-points'],  self.cache['hps-triangs']
        )

        hpsmesh.meta['edgescodes'] = self.cache['edgespaired'].copy()

        self.cache['hps-mesh'] = hpsmesh
        return hpsmesh

    def make_hps_triangs(self):

        edges = self.mesh.edges_paired()
        triangs = tables.norm_table(edges)

        self.cache['edgespaired'] = edges
        self.cache['hps-triangs'] = triangs

        return triangs

    def make_hps_points(self):

        _ = self.make_edgesnodes()
        _ = self.flat_edgesnodes()

        return _

    def make_edgesnodes(self):

        edges = 0.5 * (
            self.trinodes1 + self.trinodes2
        )

        self.cache['edges-nodes'] = edges
        return edges

    def flat_edgesnodes(self):

        points = np.zeros(
            np.amax(self.cache['hps-triangs']) + 1, dtype=complex
        )

        nodes = self.cache['edges-nodes']
        triangs = self.cache['hps-triangs']

        points[triangs.flat] = nodes.flat
        self.cache['hps-points'] = points.copy()

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
