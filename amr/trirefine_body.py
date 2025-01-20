# -*- coding: utf-8 -*-
"""Mesh body refinement.
"""
import numpy as np
from triellipt.utils import pairs, tables


def refine_mesh(mesh, trinums):
    _ = BodyRefiner.from_mesh(mesh).with_trinums_to_refine(trinums)
    return _.release_mesh()


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self.meta = self.fetch_meta()
        self.cache = {}

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    @property
    def hasvoids(self):
        return bool(
            self.meta['supvoids']
        )

    def fetch_meta(self):
        return {
            'supvoids': self.fetch_supvoids()
        }

    def fetch_supvoids(self):

        has_voids = self.mesh.hasvoids()

        if not has_voids:
            return None

        suptri = self.mesh.supertriu()
        return suptri.supvoids()


class BodyRefiner(MeshAgent):
    """Mesh body refiner.
    """

    def with_trinums_to_refine(self, trinums):
        _ = self.trinums_filter.get_filtered(trinums)
        self.meta['trinums-to-refine'] = _
        return self

    def release_mesh(self):

        mesh = self.make_mesh_beta()

        mesh = mesh.add_meta(
            {'data-refiner': self.make_data_refiner()}
        )

        return mesh

    def make_data_refiner(self):
        return self.maker_data_refiner.get_refiner()

    def make_mesh_beta(self):

        self.make_mesh_alpha()

        mesh = self.maker_mesh_beta.get_mesh()
        self.cache['mesh-beta'] = mesh
        return mesh

    def make_mesh_alpha(self):

        _ = self.make_core_mesh()

        mesh = self.maker_mesh_alpha.get_mesh()
        self.cache['mesh-alpha'] = mesh
        return mesh

    def make_core_mesh(self):
        mesh = self.maker_core_mesh.get_mesh()
        self.cache['core-mesh'] = mesh
        return mesh

    @property
    def parent_triangs(self):
        return self.mesh.triangs[self.trinums_to_refine, :]

    @property
    def parent_vertices(self):
        return self.mesh.points[self.parent_triangs]

    @property
    def trinums_to_refine(self):
        return self.meta['trinums-to-refine']

    @property
    def numbering_offset(self):
        return self.mesh.npoints

    @property
    def trinums_filter(self):
        return TrinumsFilter.from_refiner(self)

    @property
    def maker_core_mesh(self):
        return MakerCoreMesh.from_refiner(self)

    @property
    def maker_mesh_alpha(self):
        return MakerMeshAlpha.from_refiner(self)

    @property
    def maker_mesh_beta(self):
        return MakerMeshBeta.from_refiner(self)

    @property
    def maker_data_refiner(self):
        return MakerDataRefiner.from_refiner(self)


class RefinerAgent:
    """Operator on a mesh refiner.
    """

    def __init__(self, refiner):
        self.refiner = refiner
        self.meta = self.fetch_meta()
        self.cache = {}

    @classmethod
    def from_refiner(cls, refiner):
        return cls(refiner)

    def fetch_meta(self):
        return {}


class MakerCoreMesh(RefinerAgent):
    """Maker of a core mesh (pre-alpha).
    """

    def get_mesh(self):

        triangs = self.make_core_triangs()
        points = self.make_core_points()

        return self.push_core_mesh(points, triangs)

    def push_core_mesh(self, points, triangs):
        return type(self.refiner.mesh).from_data(points, triangs)

    def make_core_triangs(self):

        triangs_encoded = pairs.paircols(
            self.refiner.parent_triangs[:, [0, 1, 2, 0]]
        )

        triangs_decoded = tables.norm_table(triangs_encoded)
        core_triangs = triangs_decoded

        self.cache['core-triangs'] = core_triangs
        return core_triangs

    def make_core_points(self):

        as_table = self.make_core_points_table()
        as_array = self.push_core_points_array(as_table)

        return as_array

    def make_core_points_table(self):

        verts = self.refiner.parent_vertices

        vert1 = 0.5 * (verts[:, 0] + verts[:, 1])
        vert2 = 0.5 * (verts[:, 1] + verts[:, 2])
        vert3 = 0.5 * (verts[:, 2] + verts[:, 0])

        return _pack_cols(
            vert1, vert2, vert3
        )

    def push_core_points_array(self, points_as_table):

        points_array = np.zeros(
            self.core_points_count, dtype=complex
        )

        points_array[self.core_triangs.flat] = points_as_table.flat
        return points_array

    @property
    def core_triangs(self):
        return self.cache['core-triangs']

    @property
    def core_points_count(self):
        return np.amax(self.core_triangs) + 1


class MakerMeshAlpha(RefinerAgent):
    """Maker of an alpha-stage mesh.
    """

    EMPTY_TRIANG = np.array([[], [], []]).T

    def fetch_meta(self):
        return {
            'core-triangs-shifted': self.get_core_triangs_shifted()
        }

    def get_core_triangs_shifted(self):

        core_triangs = self.refiner.cache['core-mesh'].triangs

        return np.add(
            core_triangs, self.refiner.numbering_offset
        )

    def get_mesh(self):

        self.set_mesher_cache()
        self.run_mesher()

        return self.cache['mesh-alpha']

    def run_mesher(self):

        mesh = self.refiner.mesh

        trinums_to_del = self.cache['trinums-to-del']
        triangs_to_add = self.cache['triangs-to-add']
        points_to_add = self.cache['points-to-add']

        mesh = mesh.deltriangs(*trinums_to_del)
        mesh = mesh.add_triangs(triangs_to_add)
        mesh = mesh.add_points(points_to_add)

        self.cache['mesh-alpha'] = mesh

    def set_mesher_cache(self):
        self.cache = {
            'triangs-to-add': self.make_triangs_to_add(),
            'trinums-to-del': self.make_trinums_to_del(),
            'points-to-add': self.take_points_to_add()
        }

    def take_points_to_add(self):
        return self.refiner.cache['core-mesh'].points

    def make_trinums_to_del(self):
        return self.refiner.trinums_to_refine

    def make_triangs_to_add(self):

        cores = self.take_cores_shifted()
        sides = self.make_cores_neighbors()
        voids = self.make_voids_to_add()

        triangs = np.vstack(
            [cores, sides, voids]
        )

        return triangs.astype(int)

    def take_cores_shifted(self):
        return self.meta['core-triangs-shifted']

    def make_cores_neighbors(self):

        sup = np.hsplit(self.refiner.parent_triangs, 3)
        ker = np.hsplit(self.meta['core-triangs-shifted'], 3)

        side1 = _stack_cols(ker[0], sup[1], ker[1])
        side2 = _stack_cols(ker[1], sup[2], ker[2])
        side3 = _stack_cols(ker[2], sup[0], ker[0])

        return np.vstack(
            [side1, side2, side3]
        )

    def make_voids_to_add(self):

        pivnums, trinums, locnums = self.find_voids_pivots()

        if pivnums.size == 0:
            return self.EMPTY_TRIANG

        loc1 = (locnums + 0) % 3
        loc2 = (locnums + 1) % 3

        host_triangs = self.refiner.parent_triangs

        eastnums = host_triangs[trinums, loc1]
        westnums = host_triangs[trinums, loc2]

        voids = _pack_cols(
            eastnums, westnums, pivnums
        )

        return voids

    def find_voids_pivots(self):

        nodesmap = tables.maptable(
            self.meta['core-triangs-shifted']
        )

        nodnums, trinums, locnums = nodesmap.atrank(1)
        return nodnums, trinums, locnums


class MakerMeshBeta(RefinerAgent):
    """Maker of a beta-stage mesh.
    """

    def fetch_meta(self):
        return {
            'meshedge-alpha': self.mesh_alpha.meshedge()
        }

    def get_mesh(self):

        self.find_voids_on_edge()
        self.pull_voids_on_edge()

        return self.cache['mesh-beta']

    def find_voids_on_edge(self):

        voidsnums = self.mesh_alpha.getvoids()

        on_edge = np.isin(
            voidsnums, self.mesh_edge_alpha.trinums_unique
        )

        self.cache['voids-on-edge'] = voidsnums[on_edge]

    def pull_voids_on_edge(self):

        voids_on_edge = self.cache['voids-on-edge']

        mesh = self.mesh_alpha.deltriangs(*voids_on_edge)
        self.cache['mesh-beta'] = mesh

    @property
    def mesh_alpha(self):
        return self.refiner.cache['mesh-alpha']

    @property
    def mesh_edge_alpha(self):
        return self.meta['meshedge-alpha']


class MakerDataRefiner(RefinerAgent):
    """Makes a data-refiner for a pre-release mesh.
    """

    def get_refiner(self):
        return DataRefiner(
            self.get_nodes_images()
        )

    def get_nodes_images(self):

        old_nodes_images = self.get_old_nodes_images()
        new_nodes_images = self.get_new_nodes_images()

        return np.vstack(
            [old_nodes_images, new_nodes_images]
        )

    def get_old_nodes_images(self):

        image1 = np.arange(self.refiner.mesh.npoints)
        image2 = np.arange(self.refiner.mesh.npoints)

        return _pack_cols(image1, image2)

    def get_new_nodes_images(self):

        hosts = self.refiner.parent_triangs
        nodes = self.refiner.cache['core-mesh'].triangs

        image1 = hosts[:, [0, 1, 2]]
        image2 = hosts[:, [1, 2, 0]]

        _, ind = np.unique(nodes.flat, return_index=True)

        return _pack_cols(
            image1.flat[ind], image2.flat[ind]
        )


class TrinumsFilter(RefinerAgent):
    """Filters trinums to refine.
    """

    def get_filtered(self, trinums):

        trinums = np.array(trinums)

        trinums = self.filter_are_unique(trinums)
        trinums = self.filter_are_in_mesh(trinums)
        trinums = self.filter_out_supvoids(trinums)

        return trinums.astype(int)

    def filter_are_unique(self, trinums):
        return np.unique(trinums).astype(int)

    def filter_are_in_mesh(self, trinums):

        mask_in_mesh = np.logical_and(
            trinums >= 0, trinums < self.mesh_size
        )

        return np.compress(
            mask_in_mesh, trinums
        )

    def filter_out_supvoids(self, trinums):

        if not self.refiner.hasvoids:
            return trinums

        return np.setdiff1d(
            trinums, self.supvoids.supbodies.flat
        )

    @property
    def supvoids(self):
        return self.refiner.meta['supvoids']

    @property
    def mesh_size(self):
        return self.refiner.mesh.ntriangs


class DataRefiner:
    """Transmits data to a refined mesh.

    Attributes
    ----------
    nodes_images : two-column-int-array
        Images of refined nodes on a parent mesh.

    """

    def __init__(self, nodes_images):
        self.nodes_images = nodes_images

    def __call__(self, parent_data):
        return self.from_data_images(
            parent_data[self.nodes_images]
        )

    def from_data_images(self, data_images):
        return 0.5 * (data_images[:, 0] + data_images[:, 1])


def _stack_cols(*cols):
    return np.hstack(cols)


def _pack_cols(*cols):
    return np.vstack(cols).T.copy('C')
