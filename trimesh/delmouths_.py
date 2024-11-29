# -*- coding: utf-8 -*-
"""Removes mesh mouths.
"""
import numpy as np
from triellipt.utils import pairs


def remove_mouths(mesh):
    """Removes mesh mouths.

    Returns
    -------
    TriMesh
        Mesh without mouths.

    """

    if not mesh.hasvoids():
        return mesh

    _ = MouthsCleaner.from_mesh(mesh)
    return _.release_mesh()


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self.meta = self.fetch_meta()

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    def fetch_meta(self):

        suptri_prime = self.fetch_suptri_prime()
        suptri_voids = self.fetch_suptri_voids(suptri_prime)

        return {
            'suptri-prime': suptri_prime,
            'suptri-voids': suptri_voids
        }

    def fetch_suptri_prime(self):
        return self.mesh.supertriu()

    def fetch_suptri_voids(self, suptri_prime):
        return suptri_prime.supvoids()

    @property
    def suptri_prime(self):
        return self.meta.get('suptri-prime')

    @property
    def suptri_voids(self):
        return self.meta.get('suptri-voids')


class MouthsCleaner(MeshAgent):
    """Cleaner of mesh mouths.
    """

    def __init__(self, mesh):
        super().__init__(mesh)

        self.supmouths = None
        self.mesh_alpha = None
        self.mesh_beta = None

    def release_mesh(self):

        self.mesh_beta = self.make_mesh_beta()

        if self.mesh_beta is self.mesh:
            return self.mesh

        return self.mesh_beta.delghosts()

    def make_mesh_beta(self):

        self.mesh_alpha = self.make_mesh_alpha()

        if self.mesh_alpha is self.mesh:
            return self.mesh

        return self.maker_mesh_beta.get_mesh()

    def make_mesh_alpha(self):

        self.supmouths = self.find_supmouths()

        if self.supmouths.size == 0:
            return self.mesh

        return self.maker_mesh_alpha.get_mesh()

    def find_supmouths(self):
        return self.maker_supmouths.get_suptriu()

    @property
    def maker_supmouths(self):
        return MakerSupMouths.from_cleaner(self)

    @property
    def maker_mesh_alpha(self):
        return MakerMeshAlpha.from_cleaner(self)

    @property
    def maker_mesh_beta(self):
        return MakerMeshBeta.from_cleaner(self)


class CleanerAgent:
    """Operator on a mesh cleaner.
    """

    def __init__(self, cleaner):
        self.cleaner = cleaner

    @classmethod
    def from_cleaner(cls, cleaner):
        return cls(cleaner)


class MakerSupMouths(CleanerAgent):
    """Finds a suptriu made of mouths.
    """

    def get_suptriu(self):

        codescache = self.make_pairing_codes()
        supmouths = self.from_pairing_codes(codescache)

        return supmouths

    def make_pairing_codes(self):
        return {
            'voids_ears': self.pair_voids_ears(),
            'suptrinums': self.pair_suptrinums()
        }

    def from_pairing_codes(self, codes):

        mask = np.isin(
            codes['suptrinums'], codes['voids_ears']
        )

        inds_mouths, = np.where(
            np.sum(mask, axis=1) == 2
        )

        return self.cleaner.suptri_prime.subtriu(*inds_mouths)

    def pair_voids_ears(self):

        supv = self.cleaner.suptri_voids

        return pairs.sympaired(
            supv.trinums2, supv.trinums3
        )

    def pair_suptrinums(self):

        supt = self.cleaner.suptri_prime

        trio = _pack_cols(
            supt.trinums1, supt.trinums2, supt.trinums3
        )

        return pairs.paircols(
            trio[:, [0, 1, 2, 0]]
        )


class MakerMeshAlpha(CleanerAgent):
    """Maker of an alpha-stage mesh.
    """

    @property
    def supmouths(self):
        return self.cleaner.supmouths

    def get_mesh(self):

        meta = self.make_mesher_meta()
        mesh = self.from_mesher_meta(meta)

        return mesh

    def from_mesher_meta(self, meta):

        mesh = self.cleaner.mesh

        mesh = self.draft_new_mesh(mesh, meta)
        mesh = self.align_new_voids(mesh, meta)

        return mesh

    def align_new_voids(self, mesh, meta):

        voids = meta['triangs-to-add']['extravoids']

        west = mesh.points[voids[:, 0]]
        east = mesh.points[voids[:, 1]]

        mesh.points[voids[:, 2]] = 0.5 * (west + east)
        return mesh

    def draft_new_mesh(self, mesh, meta):

        mesh = mesh.deltriangs(*meta['trinums-to-del'])
        mesh = mesh.add_triangs(meta['triangs-to-add']['suptriangs'])
        mesh = mesh.add_triangs(meta['triangs-to-add']['extravoids'])

        return mesh

    def make_mesher_meta(self):
        return {
            'triangs-to-add': self.make_triangs_to_add(),
            'trinums-to-del': self.make_trinums_to_del()
        }

    def make_trinums_to_del(self):
        return self.supmouths.supbodies.flatten()

    def make_triangs_to_add(self):

        suptriangs = self.make_suptriangs_to_add()
        extravoids = self.make_extravoids_to_add()

        return {
            'suptriangs': suptriangs,
            'extravoids': extravoids
        }

    def make_extravoids_to_add(self):

        sup = np.hsplit(self.supmouths.supmesh.triangs, 3)
        ker = np.hsplit(self.supmouths.kermesh.triangs, 3)

        voids_0 = _stack_cols(sup[1], sup[0], ker[1])
        voids_1 = _stack_cols(sup[2], sup[1], ker[2])
        voids_2 = _stack_cols(sup[0], sup[2], ker[0])

        assert voids_0.shape[1] == 3, voids_0.shape
        assert voids_1.shape[1] == 3, voids_1.shape
        assert voids_2.shape[1] == 3, voids_2.shape

        extra_voids = np.vstack(
            [voids_0, voids_1, voids_2]
        )

        return extra_voids

    def make_suptriangs_to_add(self):
        return self.supmouths.supmesh.triangs


class MakerMeshBeta(CleanerAgent):
    """Makes a beta-mesh from an alpha-mesh.
    """

    def __init__(self, agent):
        super().__init__(agent)
        self.meta = self.fetch_meta()

    @property
    def mesh_alpha(self):
        return self.cleaner.mesh_alpha

    @property
    def voids_trinums(self):
        return self.meta.get('voids')['trinums']

    @property
    def voids_triangs(self):
        return self.meta.get('voids')['triangs']

    def fetch_meta(self):
        return {
            'voids': self.fetch_voids_meta()
        }

    def fetch_voids_meta(self):

        trinums = self.mesh_alpha.getvoids()
        triangs = self.mesh_alpha.triangs[trinums, :]

        return {
            'trinums': trinums,
            'triangs': triangs
        }

    def get_mesh(self):

        twins_nums = self.find_twins_voids()
        clean_mesh = self.pull_twins_voids(twins_nums)

        return clean_mesh

    def pull_twins_voids(self, twins_nums):

        if twins_nums.size == 0:
            return self.mesh_alpha

        return self.mesh_alpha.deltriangs(*twins_nums)

    def find_twins_voids(self):

        twins_inds = _find_twins(
            self.make_voids_codes()
        )

        return self.voids_trinums[twins_inds]

    def make_voids_codes(self):
        return pairs.sympaired(
            self.voids_triangs[:, 0], self.voids_triangs[:, 1]
        )


def _pack_cols(*cols):
    return np.vstack(cols).T.copy('C')


def _stack_cols(*cols):
    return np.hstack(cols)


def _find_twins(data):

    _, invinds, counts = np.unique(
        data, return_inverse=True, return_counts=True
    )

    twins_mask = counts[invinds] == 2

    twins_inds = np.compress(
        twins_mask, np.arange(data.size)
    )

    return twins_inds
