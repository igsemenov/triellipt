# -*- coding: utf-8 -*-
"""Public HPS unit.
"""
from triellipt.hps import hpsmesher


class HPSData:
    """Parent HPS unit.
    """

    def __init__(self, mesh=None, meta=None):
        self.mesh = mesh
        self.meta = meta

    @classmethod
    def from_mesh(cls, inputmesh):

        mesh = hpsmesher.gethpsmesh(inputmesh)

        meta = {
            'parentmesh': inputmesh,
            'edgescodes': mesh.meta['edgescodes'].copy()
        }

        mesh = mesh.twin()
        meta = meta.copy()

        return cls(mesh, meta)

    @property
    def parent_mesh(self):
        return self.meta['parentmesh']


class HPSUnit(HPSData):
    """HPS computing unit.
    """
