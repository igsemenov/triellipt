# -*- coding: utf-8 -*-
"""Mesh refinement.
"""
import numpy as np
from triellipt.utils import pairs, tables


def refine_mesh(mesh, trinums):
    _ = get_refiner(mesh, trinums)
    return _.release_mesh()


def get_refiner(mesh, trinums):
    return MeshRefiner.from_mesh(mesh).with_trinums(trinums)


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

    def fetch_meta(self):
        return {
            'voids': self.fetch_voids_meta()
        }

    def fetch_voids_meta(self):

        trinums = self.mesh.getvoids()
        triangs = self.mesh.triangs[trinums, :]

        return {
            'trinums': trinums,
            'triangs': triangs
        }

    @property
    def hasvoids(self):
        return bool(
            self.meta['voids']
        )

    @property
    def voids_trinums(self):
        return self.meta['voids']['trinums']

    @property
    def voids_triangs(self):
        return self.meta['voids']['triangs']


class MeshRefiner(MeshAgent):
    """Refines the target triangles.
    """

    def with_trinums(self, trinums):
        _ = self.trinums_filter.get_filtered(trinums)
        self.meta['trinums-to-refine'] = _
        return self

    @property
    def trinums_filter(self):
        return FilterTrinums.from_refiner(self)

    @property
    def target_trinums(self):
        return self.meta['trinums-to-refine']

    def release_mesh(self):

        mesh = self.make_mesh_beta()

        mesh = mesh.twin()
        mesh = mesh.add_meta(
            {'data-refiner': self.make_data_refiner()}
        )

        return mesh

    def make_data_refiner(self):
        return self.maker_data_refiner.get_refiner()

    def make_mesh_beta(self):

        _ = self.make_mesh_alpha()
        _ = self.maker_mesh_beta.get_mesh()

        self.cache['mesh-beta'] = _
        return _

    def make_mesh_alpha(self):

        _ = self.make_core_mesh()
        _ = self.maker_mesh_alpha.get_mesh()

        self.cache['mesh-alpha'] = _
        return _

    def make_core_mesh(self):
        _ = self.maker_core_mesh.get_mesh()
        self.cache['core-mesh'] = _
        return _

    @property
    def parent_triangs(self):
        return self.mesh.triangs[self.target_trinums, :]

    @property
    def parent_vertices(self):
        return self.mesh.points[self.parent_triangs]

    @property
    def numbering_offset(self):
        return self.mesh.npoints

    @property
    def maker_core_mesh(self):
        return MakerMeshCore.from_refiner(self)

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

    @property
    def target_submesh(self):
        return self.refiner.target_submesh


class MakerMeshCore(RefinerAgent):
    """Mesh from midpoints of target-triangles-edges.
    """

    def get_mesh(self):

        triangs = self.make_triangs()
        points = self.make_points()

        return self.push_core_mesh(points, triangs)

    def push_core_mesh(self, points, triangs):
        return self.refiner.mesh.from_data(points, triangs)

    def make_triangs(self):

        triangs_encoded = self.target_submesh.edges_paired()
        triangs_decoded = tables.norm_table(triangs_encoded)

        core_triangs = triangs_decoded
        self.cache['core-triangs'] = core_triangs

        return core_triangs

    def make_points(self):

        as_table = self.make_points_table()
        as_array = self.push_points_array(as_table)

        return as_array

    def make_points_table(self):

        verts = self.refiner.parent_vertices

        vert1 = 0.5 * (verts[:, 0] + verts[:, 1])
        vert2 = 0.5 * (verts[:, 1] + verts[:, 2])
        vert3 = 0.5 * (verts[:, 2] + verts[:, 0])

        return _pack_cols(
            vert1, vert2, vert3
        )

    def push_points_array(self, points_as_table):

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

    @property
    def target_submesh(self):
        return self.refiner.mesh.submesh(*self.target_trinums)

    @property
    def target_trinums(self):
        return self.refiner.target_trinums


class MakerMeshAlpha(RefinerAgent):
    """Extends the core mesh.

    - Adds neighbours around core triangles.
    - Adds voids for one-rank new points.

    """

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
            'points-to-add': self.make_points_to_add()
        }

    def make_points_to_add(self):
        return self.refiner.cache['core-mesh'].points

    def make_trinums_to_del(self):
        return self.refiner.target_trinums

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

        sup = np.hsplit(
            self.refiner.parent_triangs, 3
        )

        ker = np.hsplit(
            self.meta['core-triangs-shifted'], 3
        )

        side1 = _stack_cols(ker[0], sup[1], ker[1])
        side2 = _stack_cols(ker[1], sup[2], ker[2])
        side3 = _stack_cols(ker[2], sup[0], ker[0])

        return np.vstack(
            [side1, side2, side3]
        )

    def make_voids_to_add(self):

        pivnums, trinums, locnums = self.find_new_voids_pivots()

        if pivnums.size == 0:
            return np.array([[], [], []]).T

        loc1 = (locnums + 0) % 3
        loc2 = (locnums + 1) % 3

        host_triangs = self.refiner.parent_triangs

        eastnums = host_triangs[trinums, loc1]
        westnums = host_triangs[trinums, loc2]

        voids = _pack_cols(
            eastnums, westnums, pivnums
        )

        return voids

    def find_new_voids_pivots(self):

        nodesmap = tables.maptable(
            self.meta['core-triangs-shifted']
        )

        nodnums, trinums, locnums = nodesmap.atrank(1)
        return nodnums, trinums, locnums


class MakerMeshBeta(RefinerAgent):
    """Cleans the alpha-mesh.

    - Removes twin voids.
    - Removes voids from the edge.

    """

    def fetch_meta(self):
        return {
            'mesh-alpha-edge': self.mesh_alpha.meshedge(),
            'mesh-alpha-voids': self.mesh_alpha_voids_meta()
        }

    def mesh_alpha_voids_meta(self):

        trinums = self.mesh_alpha.getvoids()
        triangs = self.mesh_alpha.triangs[trinums, :]

        return {
            'trinums': trinums,
            'triangs': triangs
        }

    @property
    def mesh_alpha(self):
        return self.refiner.cache['mesh-alpha']

    def get_mesh(self):
        """Creates the beta-mesh.
        """

        _ = self.set_mesher_cache()
        _ = self.run_mesher()

        return _

    def run_mesher(self):

        self.cache['mesh-beta'] = self.mesh_alpha

        _ = self.remove_voids()
        _ = self.reconnect_twins_pivs()

        return _

    def remove_voids(self):

        mesh = self.cache['mesh-beta']

        mesh = mesh.deltriangs(
            *self.get_voids_to_remove()
        )

        self.cache['mesh-beta'] = mesh
        return mesh

    def get_voids_to_remove(self):

        voids_on_edge = self.cache['voids-on-edge']

        if 'voids-twins' not in self.cache:
            return voids_on_edge

        return np.r_[
            voids_on_edge, self.cache['voids-twins']['nums'].flat
        ]

    def reconnect_twins_pivs(self):

        mesh = self.cache['mesh-beta']

        if 'voids-twins' not in self.cache:
            return mesh

        pivs1 = self.cache['voids-twins']['pivs'][:, 0]
        pivs2 = self.cache['voids-twins']['pivs'][:, 1]

        mesh = _reconnect_nodes(mesh, pivs1, pivs2)
        mesh = mesh.delghosts()

        meta = {
            'nodes-images-alpha': mesh.meta['old-nodes-numbers']
        }

        mesh = mesh.twin()
        mesh = mesh.add_meta(meta)

        return mesh

    def set_mesher_cache(self):

        _ = self.find_voids_on_edge()
        _ = self.find_voids_twins()

        return True

    def find_voids_twins(self):

        trinums = self.meta['mesh-alpha-voids']['trinums']
        triangs = self.meta['mesh-alpha-voids']['triangs']

        voids_codes = pairs.sympaired(
            triangs[:, [0]], triangs[:, [1]]
        )

        voids_map = tables.maptable(voids_codes)
        _, twin_inds, _ = voids_map.atrank(2)

        if twin_inds.size == 0:
            return {}

        twin_inds = _as2col(twin_inds)
        twin_nums = trinums[twin_inds]

        twin_pivs = _stack_cols(
            triangs[twin_inds[:, [0]], 2],
            triangs[twin_inds[:, [1]], 2]
        )

        twin_pivs = np.sort(
            twin_pivs, axis=1
        )

        twin_meta = {
            'nums': twin_nums,
            'pivs': twin_pivs
        }

        self.cache['voids-twins'] = twin_meta
        return twin_meta

    def find_voids_on_edge(self):

        edge = self.meta['mesh-alpha-edge']

        trinums = self.meta['mesh-alpha-voids']['trinums']
        triangs = self.meta['mesh-alpha-voids']['triangs']

        voids_codes = pairs.sympaired(
            triangs[:, 0], triangs[:, 1]
        )

        edge_codes = pairs.sympaired(
            edge.nodnums1, edge.nodnums2
        )

        mask = np.isin(
            voids_codes, edge_codes
        )

        where_mask, = np.where(mask)
        voids_on_edge = trinums[where_mask]

        self.cache['voids-on-edge'] = voids_on_edge
        return voids_on_edge


class MakerDataRefiner(RefinerAgent):
    """Makes a data-refiner for a pre-release mesh.
    """

    def get_refiner(self):

        meta = {
            'nodes-images': self.get_nodes_images_beta()
        }

        refiner = DataRefiner(
            self.refiner.mesh, meta
        )

        return refiner

    def get_nodes_images_beta(self):
        """Nodes images for the beta-mesh.
        """

        images_alpha = self.get_nodes_images_alpha()
        images_beta_in_alpha = self.nodes_images_beta_in_alpha()

        if images_beta_in_alpha is None:
            return images_alpha.copy('C')

        images_beta = images_alpha[images_beta_in_alpha, :]
        return images_beta.copy('C')

    def nodes_images_beta_in_alpha(self):
        beta = self.refiner.cache['mesh-beta']
        return beta.meta.get('nodes-images-alpha')

    def get_nodes_images_alpha(self):
        """Nodes images for the alpha-mesh.
        """

        old_nodes_images = self.get_old_nodes_images()
        new_nodes_images = self.get_new_nodes_images()

        images_alpha = np.vstack(
            [old_nodes_images, new_nodes_images]
        )

        return images_alpha

    def get_old_nodes_images(self):

        image1 = np.arange(self.refiner.mesh.npoints)
        image2 = np.arange(self.refiner.mesh.npoints)

        return _pack_cols(image1, image2)

    def get_new_nodes_images(self):

        hosts = self.refiner.parent_triangs
        nodes = self.refiner.cache['core-mesh'].triangs

        image1 = hosts[:, [0, 1, 2]]
        image2 = hosts[:, [1, 2, 0]]

        _, ind = np.unique(
            nodes.flat, return_index=True
        )

        return _pack_cols(
            image1.flat[ind], image2.flat[ind]
        )


class FilterTrinums(RefinerAgent):
    """Filters trinums to refine.
    """

    def get_filtered(self, trinums):

        trinums = np.array(trinums)

        trinums = self.filter_are_unique(trinums)
        trinums = self.filter_are_in_mesh(trinums)
        trinums = self.filter_touching_voids(trinums)

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

    def filter_touching_voids(self, trinums):

        if not self.refiner.hasvoids:
            return trinums

        mask = np.isin(
            self.codes_target_submesh(trinums), self.codes_voids_ears()
        )

        return np.delete(
            trinums, np.where(np.any(mask, axis=1))[0]
        )

    def codes_voids_ears(self):

        voids_triangs = self.refiner.voids_triangs

        west_ears_codes = pairs.sympaired(
            voids_triangs[:, 1], voids_triangs[:, 2]
        )

        east_ears_codes = pairs.sympaired(
            voids_triangs[:, 0], voids_triangs[:, 2]
        )

        return np.hstack(
            [west_ears_codes, east_ears_codes]
        )

    def codes_target_submesh(self, trinums):
        _ = self.refiner.mesh.submesh(*trinums)
        return _.edges_paired()

    @property
    def mesh_size(self):
        return self.refiner.mesh.ntriangs


class DataRefiner:
    """Transmits data to a refined mesh.

    Attributes
    ----------
    mesh : TriMesh
        Source triangle mesh.
    meta : dict
        Refiner meta-data.

    """

    def __init__(self, mesh, meta):
        self.mesh = mesh
        self.meta = meta

    @property
    def nodes_images(self):
        """Images of the refined nodes in the source mesh.
        """
        return self.meta['nodes-images']

    def __call__(self, data):
        return self.refine(data)

    def refine(self, source_data):
        return self.from_data_images(
            source_data[self.nodes_images]
        )

    def from_data_images(self, data_images):
        return 0.5 * (
            data_images[:, 0] + data_images[:, 1]
        )

    @property
    def master_nodes(self):
        return self.nodes_range[
            self.nodes_images[:, 0] == self.nodes_images[:, 1]
        ]

    @property
    def slaves_nodes(self):
        return self.nodes_range[
            self.nodes_images[:, 0] != self.nodes_images[:, 1]
        ]

    @property
    def nodes_range(self):
        return np.arange(
            self.nodes_images.shape[0]
        )


def _stack_cols(*cols):
    return np.hstack(cols)


def _pack_cols(*cols):
    return np.vstack(cols).T.copy('C')


def _as2col(data):
    return data.reshape(data.size // 2, 2)


def _reconnect_nodes(mesh, nodes1, nodes2):
    """Replaces the nodes2 by the nodes1 in a mesh.
    """

    newnodes = np.arange(mesh.npoints)
    newnodes[nodes2] = nodes1

    new_triangs = newnodes[mesh.triangs]
    return mesh.update_triangs(new_triangs)


def constr_data(mesh, data):
    """Constrains data on a mesh.
    """

    voids_trinums = mesh.getvoids()

    if voids_trinums.size == 0:
        return data
    voids_triangs = mesh.triangs[voids_trinums, :]

    west, east, pivs = voids_triangs.T

    data[pivs] = 0.5 * (
        data[west] + data[east]
    )

    return data.copy('C')
