# -*- coding: utf-8 -*-
"""Global mesh skeleton.
"""
import numpy as np
from triellipt.fem.femoprs import mesh_geom
from triellipt.mrs import premesh


def getskel(unit):
    """Makes the MRS skeleton.
    """
    return TriSkel.from_mrs_unit(unit)


class TriSkel:
    """Skeleton of the global mesh.
    """

    def __init__(self):
        self.unit = None
        self.body = None
        self.cache = {}

    @classmethod
    def from_mrs_unit(cls, unit):
        return TriSkelMaker.from_mrs_unit(unit).get_skel()

    @property
    def data(self):
        return self.body['data']

    @property
    def areas(self):
        if 'areas' in self.cache:
            return self.cache['areas']
        self.cache['areas'] = self.get_host_areas()
        return self.areas

    @property
    def ownmesh(self):
        if 'ownmesh' in self.cache:
            return self.cache['ownmesh']
        self.cache['ownmesh'] = self.asmesh()
        return self.ownmesh

    @property
    def hostmesh(self):
        return self.unit.mesh

    @property
    def triangs(self):
        return self.body['triangs']

    @property
    def triedges(self):
        return self.body['triedges']

    @property
    def ranks(self):
        return self.unit.ranks

    def asmesh(self):
        return TriSkelMesher(self).get_mesh()

    def get_host_areas(self):
        return mesh_geom(self.hostmesh).areas

    def get_premesh(self):
        return premesh.getpremesh(self)

    def get_premesher(self):
        return premesh.getpremesher(self)


class MRSUnitAgent:
    """Operator on an MRS unit.
    """

    def __init__(self, unit):
        self.unit = unit
        self.meta = self.fetch_meta()

    @classmethod
    def from_mrs_unit(cls, unit):
        return cls(unit)

    def fetch_meta(self):

        numshifts = np.add.accumulate(
            np.r_[0, self.unit.counts[:-1]]
        )

        return {
            'numshifts': numshifts
        }

    @property
    def shifts(self):
        return self.meta['numshifts']


class TriSkelMaker(MRSUnitAgent):
    """Makes the global mesh skeleton.
    """

    def get_skel(self):

        body = self.make_body()
        skel = self.from_body(body)

        return skel

    def from_body(self, body):

        skel = TriSkel()

        skel.unit = self.unit
        skel.body = body

        return skel

    def make_body(self):
        return {
            'data': self.stack_data(),
            'points': self.stack_points(),
            'triangs': self.stack_triangs(),
            'triedges': self.stack_triedges()
        }

    def stack_triangs(self):

        zip_elms_shifts = zip(
            self.unit.elms, self.shifts
        )

        chunks = [
            e.mesh.triangs + s for e, s in zip_elms_shifts
        ]

        return np.vstack(chunks)

    def stack_points(self):

        chunks = [
            e.mesh.points for e in self.unit.elms
        ]

        return np.hstack(chunks)

    def stack_data(self):

        if not self.unit.data_keys:
            return {}

        data = {
            k: self.stack_data_at_key(k) for k in self.unit.data_keys
        }

        return data

    def stack_data_at_key(self, key):

        chunks = [
            e.data[key] for e in self.unit.elms
        ]

        return np.hstack(chunks)

    def stack_triedges(self):

        zip_elms_shifts = zip(
            self.unit.elms, self.shifts
        )

        data = [
            e.get_edges(s) for e, s in zip_elms_shifts
        ]

        return TriPacks(
            np.hstack(data), self.unit.ranks
        )


class TriSkelMesher:
    """Makes a mesh from the MRS skel.
    """

    def __init__(self, skel):
        self.skel = skel

    def get_mesh(self):

        data = self.take_data()
        points = self.take_points()
        triangs = self.take_triangs()

        mesh = self.skel.unit.mesh.from_data(points, triangs)

        mesh.meta = {
            'data': data.copy()
        }

        return mesh

    def take_data(self):
        return self.skel.body['data']

    def take_triangs(self):
        return self.skel.body['triangs']

    def take_points(self):
        return self.skel.body['points']


class TriPacks:
    """Power-of-two-sized triplets.
    """

    def __init__(self, data, ranks):
        self.data = data
        self.ranks = ranks
        self.meta = self.fetch_meta()

    @property
    def shifts(self):
        return self.meta['shifts']

    @property
    def chunks(self):
        return self.meta['chunks']

    def fetch_meta(self):
        return {
            'chunks': self.get_chunklens(),
            'shifts': self.get_numshifts()
        }

    def get_numshifts(self):

        pack_sizes = 3 * self.get_chunklens()

        return np.add.accumulate(
            np.r_[0, pack_sizes[:-1]]
        )

    def get_chunklens(self):
        return 1 + (2 ** self.ranks)

    def __call__(self, trinum, locnum):

        shift = self.shifts[trinum]
        chunk = self.chunks[trinum]

        pos1 = shift + locnum * chunk
        pos2 = shift + locnum * chunk + chunk

        return self.data[pos1:pos2]
