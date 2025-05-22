# -*- coding: utf-8 -*-
"""Global MRS mesh.
"""
import numpy as np


def getmesher(premesh):
    return MRSMesher.from_premesh(premesh)


class MRSMesh:
    """Global MRS mesh.
    """

    def __init__(self):
        self.mesh = None
        self.data = None
        self.meta = None

    @classmethod
    def from_premesh(cls, premesh):
        return ...


class PreMeshAgent:
    """Operator on a premesh.
    """

    def __init__(self, premesh):
        self.prem = premesh
        self.meta = self.fetch_meta()

    @classmethod
    def from_premesh(cls, premesh):
        return cls(premesh)

    def fetch_meta(self):
        return {
            'areas': self.fetch_areas(),
            'ranks': self.fetch_ranks()
        }

    def fetch_areas(self):

        dat = self.prem.skel.areas
        pad = self.pad_spec_triangs

        return np.pad(
            dat, pad_width=pad, mode='constant', constant_values=0
        )

    def fetch_ranks(self):

        dat = self.prem.skel.ranks
        pad = self.pad_spec_triangs

        return np.pad(
            dat, pad_width=pad, mode='constant', constant_values=1
        )

    @property
    def areas(self):
        return self.meta['areas']

    @property
    def ranks(self):
        return self.meta['ranks']

    @property
    def pad_spec_triangs(self):
        return (
            (0, self.pad_width_triangs)
        )

    @property
    def pad_width_triangs(self):
        return self.prem.mesh.size - self.prem.skel.ranks.size


class MRSMesher(PreMeshAgent):
    """Creates the MRS mesh.
    """

    def get_mrsmesh(self):
        pass
