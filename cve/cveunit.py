# -*- coding: utf-8 -*-
"""Public CVE unit.
"""
from triellipt.cve import cvemesher


def getunit(mesh):
    return CVEUnit.from_mesh(mesh)


class CVEData:
    """Parent CVE unit.
    """

    def __init__(self, mesh=None, meta=None):
        self.mesh = mesh
        self.meta = meta

    @classmethod
    def from_mesh(cls, inputmesh):

        mesh = cvemesher.getcvemesh(inputmesh)

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


class CVEUnit(CVEData):
    """CVE computing unit.
    """
