# -*- coding: utf-8 -*-
"""Removes mesh mouths.
"""
import numpy as np


def clean_mesh(mesh):
    """Removes mesh mouths.

    Returns
    -------
    TriMesh
        Mesh without mouths.

    """

    if not mesh.hasvoids():
        return mesh.twin()

    _ = make_cleaner(mesh)
    return _.mesh_cleaned()


def make_cleaner(mesh):

    cleaner = DelMouthes(mesh)

    cleaner = cleaner.enriched()
    cleaner = cleaner.with_tongues()
    cleaner = cleaner.with_mouthes()

    return cleaner


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh=None):
        self.mesh = mesh
        self.is_enriched = False

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    @property
    def mesh_twin(self):
        return self.mesh.twin()

    def enriched(self):

        self.add_voids_nums()
        self.add_super_triu()

        self.is_enriched = True
        return self

    def add_voids_nums(self):
        self.mesh.voids_nums = self.mesh.getvoids()

    def add_super_triu(self):
        self.mesh.super_triu = self.mesh.supertriu()

    @property
    def voids_trinums(self):
        return self.mesh.voids_nums

    @property
    def voids_triangs(self):
        return self.mesh.triangs[self.voids_trinums, :]

    @property
    def voids_pivots(self):
        return self.mesh.triangs[self.voids_trinums, 2]


class GetMouthes(MeshAgent):
    """Finds mesh mouths.
    """

    def __init__(self, mesh=None):
        super().__init__(mesh)

        self.tongues = None
        self.mouthes = None

    def with_tongues(self):
        self.tongues = self.find_tongues()
        return self

    def with_mouthes(self):
        self.mouthes = self.make_mouthes()
        return self

    @property
    def is_empty(self):
        if self.tongues and self.mouthes:
            return False
        return True

    def find_tongues(self):

        if self.is_enriched is False:
            return None

        return GetTongues(self.mesh).find_tongues()

    def make_mouthes(self):

        if not self.tongues:
            return None

        return self.mesh.super_triu.atcores(*self.tongues.trinums)

    def find_triangs_to_del(self):

        voids = self.pick_voids_to_del()
        cores = self.pick_cores_to_del()

        return np.hstack(
            [voids, cores]
        )

    def pick_voids_to_del(self):

        mask = np.isin(
            self.voids_pivots, self.tongues.triangs
        )

        return self.voids_trinums[mask]

    def pick_cores_to_del(self):
        return self.mouthes.supbodies.flatten()

    def make_triangs_to_add(self):

        new_cores = self.pick_new_cores()
        new_voids = self.make_new_voids()

        return {
            'new-cores': new_cores.astype(int),
            'new-voids': new_voids.astype(int)
        }

    def pick_new_cores(self):
        return self.mouthes.supmesh.triangs

    def make_new_voids(self):

        vtx0, vtx1 = self.pick_new_voids_bases()
        vtx2 = self.pick_new_voids_pivots()

        triangs = _pack_cols(vtx0, vtx1, vtx2)
        return triangs

    def pick_new_voids_bases(self):

        suptris = self.mouthes.supmesh.triangs
        locnums = self.tongues.locnums_apexes

        locnums1 = (locnums + 0) % 3
        locnums2 = (locnums + 2) % 3

        yield _pick_one_per_row(suptris, locnums1)
        yield _pick_one_per_row(suptris, locnums2)

    def pick_new_voids_pivots(self):
        return self.tongues.nodnums_apexes


class DelMouthes(GetMouthes):
    """Removes mesh mouthes.
    """

    def mesh_cleaned(self):

        if self.is_empty:
            return self.mesh_twin

        tridata = self.pick_triangs_data()
        newmesh = self.from_triangs_data(tridata)

        return newmesh

    def pick_triangs_data(self):
        return {
            'trinums-to-del': self.find_triangs_to_del(),
            'triangs-to-add': self.make_triangs_to_add()
        }

    def from_triangs_data(self, tridata):

        mesh = self.mesh

        mesh = self.del_mesh_triangs(
            mesh, tridata['trinums-to-del']
        )

        mesh = self.add_mesh_triangs(
            mesh, tridata['triangs-to-add']
        )

        mesh = self.del_mesh_ghosts(mesh)
        return mesh

    def del_mesh_ghosts(self, mesh):
        return mesh.delghosts()

    def del_mesh_triangs(self, mesh, trinums_to_del):
        return mesh.deltriangs(*trinums_to_del)

    def add_mesh_triangs(self, mesh, triangs_to_add):

        new_cores = triangs_to_add['new-cores']
        new_voids = triangs_to_add['new-voids']

        new_triangs = np.vstack(
            [new_cores, new_voids]
        )

        mesh = mesh.add_triangs(new_triangs)

        mesh = self.flat_new_voids(mesh, new_voids)
        return mesh

    def flat_new_voids(self, mesh, voids_triangs):

        new_points = _flat_triplets(
            mesh.points, voids_triangs.T
        )

        return mesh.update_points(new_points)


class GetTongues(MeshAgent):
    """Finds tongues of mesh mouths.
    """

    def find_tongues(self):
        """Returns toungues (hangers that are not corners).
        """

        hangers, corners = self.find_hangers_corners()
        tongues = self.diff_hangers_corners(hangers, corners)

        if tongues.size == 0:
            return None
        return tongues

    def diff_hangers_corners(self, hangers, corners):
        return hangers.deltriangs(corners)

    def find_hangers_corners(self):
        yield self.find_hangers()
        yield self.find_corners()

    def find_hangers(self):
        return GetHangers(self.mesh).find_hangers()

    def find_corners(self):
        return GetCorners(self.mesh).find_corners()


class GetHangers(MeshAgent):
    """Finds triangles with two voids pivots as vertices.
    """

    def find_hangers(self):

        mask = self.mask_triangs_in_pivots()
        data = self.from_triangs_in_pivots(mask)

        return data

    def from_triangs_in_pivots(self, mask):

        rows, cols = _pick_rows_one_true(mask)

        trinums_two_pivots = rows
        locnums_not_pivots = cols

        data = np.vstack(
            [trinums_two_pivots, locnums_not_pivots]
        )

        return Hangers(
            self.mesh, data
        )

    def mask_triangs_in_pivots(self):
        return np.isin(
            self.mesh.triangs, self.voids_pivots, invert=True
        )


class GetCorners(MeshAgent):
    """Finds triangles with two voids as neighbours.
    """

    @property
    def suptri_bodies(self):
        return np.vstack(
            list(self.gen_suptri_bodies())
        )

    def gen_suptri_bodies(self):
        yield self.mesh.super_triu.trinums1
        yield self.mesh.super_triu.trinums2
        yield self.mesh.super_triu.trinums3

    def find_corners(self):
        """Returns numbers of triangles that are corners.
        """

        mask = self.mask_trinums_touch_two_voids()
        nums = self.pick_trinums_touch_two_voids(mask)

        return nums

    def pick_trinums_touch_two_voids(self, mask):
        return self.mesh.super_triu.trinums[mask]

    def mask_trinums_touch_two_voids(self):
        return _mask_cols_two_trues(
            self.mask_suptriangs_in_voids()
        )

    def mask_suptriangs_in_voids(self):
        return np.isin(
            self.suptri_bodies, self.voids_trinums
        )


class Hangers:
    """Triangles with two voids-pivots as vertices. 
    """

    def __init__(self, mesh=None, data=None):
        self.mesh = mesh
        self.data = data

    @property
    def size(self):
        return self.data.shape[1]

    def update_data(self, newdata):
        return self.__class__(
            self.mesh, newdata.copy('C')
        )

    def deltriangs(self, *trinums):
        """Deletes hangers with the specified global numbers.
        """

        mask = np.isin(
            self.trinums, trinums, invert=True
        )

        return self.update_data(
            newdata=self.data[:, mask]
        )

    @property
    def trinums(self):
        return self.data[0, :]

    @property
    def triangs(self):
        """Triangles table for hangers.
        """
        return self.mesh.triangs[self.trinums, :]

    @property
    def locnums_apexes(self):
        """Local numbers of non-pivot vertices in hangers. 
        """
        return self.data[1, :]

    @property
    def nodnums_apexes(self):
        return self.mesh.triangs[
            self.trinums, self.locnums_apexes
        ]


def _mask_cols_two_trues(mask_in_rows):
    return np.sum(mask_in_rows, axis=0) == 2


def _pick_rows_one_true(mask_table):

    rows, cols = np.where(mask_table)

    _, retinds, counts = np.unique(
        rows, return_index=True, return_counts=True
    )

    one_true = retinds[counts == 1]

    yield rows[one_true]
    yield cols[one_true]


def _pick_one_per_row(table, colinds):
    return table[
        np.arange(table.shape[0]), colinds
    ]


def _pack_cols(*cols):
    return np.vstack(list(cols)).T.copy('C')


def _flat_triplets(points, triplets):

    pos0, pos1, pos2 = triplets

    points[pos2] = 0.5 * (
        points[pos0] + points[pos1]
    )

    return points
