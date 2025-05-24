# -*- coding: utf-8 -*-
"""Global MRS mesh.
"""
import numpy as np


def getmrsmesh(unit):
    return MRSMesh.from_premesh(unit.getpremesh())


class MRSMesh:
    """Global MRS mesh.
    """

    def __init__(self, prem=None, mesh=None, data=None):
        self.prem = prem
        self.mesh = mesh
        self.data = data

    @classmethod
    def from_premesh(cls, prem):

        mesh = MRSMeshMaker.from_premesh(prem).get_mrsmesh()

        return cls(
            prem, mesh, {}
        )


class PreMeshAgent:
    """Operator on a pre-mesh.
    """

    def __init__(self, premesh):
        self.prem = premesh

    @classmethod
    def from_premesh(cls, premesh):
        return cls(premesh)


class MRSMeshMaker(PreMeshAgent):
    """Makes the mesh from premesh.
    """

    def get_mrsmesh(self):

        mesh = self.prem.mesh.twin()

        mesh = self.merge_holes(mesh)
        mesh = self.merge_gluers(mesh)

        mesh = mesh.delghosts()
        inds = mesh.meta['old-nodes-numbers'].copy('C')

        mesh = mesh.twin()

        mesh.meta = {
            'nodes-in-premesh': inds
        }

        return mesh

    def merge_holes(self, mesh):

        for hole in self.prem.holes:

            pivot = hole[0]
            stack = hole[1::]

            mask = np.isin(
                mesh.triangs, stack
            )

            mesh.triangs[mask] = pivot

        return mesh

    def merge_gluers(self, mesh):
        return mesh
