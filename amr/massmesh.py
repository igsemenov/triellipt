# -*- coding: utf-8 -*-
"""Mass-collector mesh.
"""
import numpy as np
from triellipt.utils import pairs, tables


def get_massmesh(suptri):
    _ = get_mesher(suptri)
    return _.release_mesh()


def get_mesher(suptri):
    return MassMesher.from_suptri(suptri)


def get_triangs_alpha(suptri):
    _ = get_mesher(suptri)
    return _.make_triangs_alpha()


def get_triangs_beta(suptri):
    _ = get_mesher(suptri)
    return _.make_triangs_beta()


def get_mesh_gamma(suptri):
    _ = get_mesher(suptri)
    return _.make_mesh_gamma()


class SupTriAgent:
    """Operator on a supertriu.
    """

    def __init__(self, suptri):
        self.suptri = suptri
        self.meta = self.fetch_meta()
        self.cache = {}

    def fetch_meta(self):
        return {
            'root-mesh-edge': self.fetch_mesh_edge(),
            'root-mesh-voids': self.fetch_voids_meta()
        }

    def fetch_mesh_edge(self):
        return self.suptri.mesh.meshedge()

    def fetch_voids_meta(self):

        trinums = self.root_mesh.getvoids()
        triangs = self.root_mesh.triangs[trinums, :]

        return {
            'trinums': trinums, 'pivots': triangs[:, 2]
        }

    @property
    def root_mesh(self):
        return self.suptri.mesh

    @classmethod
    def from_suptri(cls, suptri):
        return cls(suptri)


class MassMesher(SupTriAgent):
    """Maker of the mass-mesh.
    """

    def release_mesh(self):

        _ = self.make_mesh_gamma()
        _ = self.make_mesh_release()

        return _

    def make_mesh_release(self):
        return self.maker_mesh_release.get_mesh()

    def make_mesh_gamma(self):

        _ = self.make_triangs_beta()

        _ = self.maker_mesh_gamma.get_mesh()
        self.cache['mesh-gamma'] = _
        return _

    def make_triangs_beta(self):

        _ = self.make_triangs_alpha()

        _ = self.maker_triangs_beta.get_triangs()
        self.cache['triangs-beta'] = _
        return _

    def make_triangs_alpha(self):
        _ = self.maker_triangs_alpha.get_table()
        self.cache['triangs-alpha'] = _
        return _

    @property
    def maker_triangs_alpha(self):
        return MakerTriangsAlpha.from_agent(self)

    @property
    def maker_triangs_beta(self):
        return MakerTriangsBeta.from_agent(self)

    @property
    def maker_mesh_gamma(self):
        return MakerMeshGamma.from_agent(self)

    @property
    def maker_mesh_release(self):
        return MakerMeshRelease.from_agent(self)


class SubSupTriAgent:
    """Operator on a supertriu agent.
    """

    def __init__(self, agent):
        self.agent = agent
        self.meta = self.fetch_meta()
        self.cache = {}

    def fetch_meta(self):
        return {}

    @classmethod
    def from_agent(cls, agent):
        return cls(agent)


class MakerTriangsAlpha(SubSupTriAgent):
    """Makes the alpha-triangs from the input supertriu.

    - Makes inner voids with the original numbers.

    """

    def fetch_meta(self):
        return {
            'sup-vertices': self.fetch_sup_vertices(),
            'ker-vertices': self.fetch_ker_vertices()
        }

    def fetch_sup_vertices(self):
        return np.hsplit(
            self.agent.suptri.supmesh.triangs, 3
        )

    def fetch_ker_vertices(self):
        return np.hsplit(
            self.agent.suptri.kermesh.triangs, 3
        )

    def get_table(self):
        """Creates the triangs-alpha-table.
        """

        table = np.copy(
            self.make_voids_vstack(), order='C'
        )

        return table.astype(int)

    def make_voids_vstack(self):

        hvoids = self.make_voids_hstack()

        return hvoids.reshape(
            (3 * self.agent.suptri.size, 3)
        )

    def make_voids_hstack(self):
        return _stack_cols(
            *self.make_voids()
        )

    def make_voids(self):

        sup = self.meta['sup-vertices']
        ker = self.meta['ker-vertices']

        voids0 = _stack_cols(sup[0], sup[1], ker[1])
        voids1 = _stack_cols(sup[1], sup[2], ker[2])
        voids2 = _stack_cols(sup[2], sup[0], ker[0])

        return voids0, voids1, voids2


class MakerTriangsBeta(SubSupTriAgent):
    """Makes the beta-triangs from the alpha-triangs.

    - Makes inner voids with the shifted numbers.
    - Makes inner cores with the shifted numbers.
    - Makes inner sides with the shifted numbers.

    - Makes extra data on voids:
        - Makes outer voids to insert.
        - Finds original voids to delete.

    """

    def fetch_meta(self):
        return {
            'triangs-alpha-codes': self.fetch_triangs_alpha_codes(),
        }

    def fetch_triangs_alpha_codes(self):
        return pairs.sympaired(
            self.triangs_alpha[:, [0]], self.triangs_alpha[:, [1]]
        )

    @property
    def triangs_alpha(self):
        return self.agent.cache['triangs-alpha']

    @property
    def triangs_alpha_size(self):
        return self.triangs_alpha.shape[0]

    @property
    def numbers_offset(self):
        return self.agent.suptri.mesh.npoints

    def get_triangs(self):
        return {
            'inner-voids': self.make_inner_voids(),
            'inner-cores': self.make_inner_cores(),
            'inner-sides': self.make_inner_sides(),
            'extra-voids': self.make_extra_voids_meta()
        }

    def make_inner_voids(self):

        triangs = np.copy(
            self.triangs_alpha, order='C'
        )

        triangs[:, 2] = self.make_inner_pivots()

        return np.copy(
            triangs, order='C'
        )

    def make_inner_cores(self):

        pivs = self.cache['inner-pivots']

        cores = np.vstack(
            np.split(pivs, pivs.size // 3)
        )

        self.cache['inner-cores'] = cores
        return cores

    def make_inner_sides(self):

        ker = np.hsplit(
            self.cache['inner-cores'], 3
        )

        sup = np.hsplit(
            self.agent.suptri.supmesh.triangs, 3
        )

        side1 = _stack_cols(sup[0], ker[0], ker[2])
        side2 = _stack_cols(sup[1], ker[1], ker[0])
        side3 = _stack_cols(sup[2], ker[2], ker[1])

        return np.vstack(
            [side1, side2, side3]
        )

    def make_inner_pivots(self):

        pivots = np.arange(
            self.triangs_alpha_size
        )

        pivots = pivots + self.numbers_offset

        self.cache['inner-pivots'] = pivots
        return pivots

    def make_extra_voids_meta(self):

        _ = self.find_alpha_voids_hanging()

        return {
            'trinums-to-del': self.find_root_voids_to_del(),
            'triangs-to-add': self.make_extra_voids_to_add()
        }

    def find_root_voids_to_del(self):

        root_voids = self.agent.meta['root-mesh-voids']
        alpha_voids = self.cache['alpha-voids-hanging']

        double_voids = np.isin(
            root_voids['pivots'], alpha_voids['pivots']
        )

        return root_voids['trinums'][double_voids]

    def make_extra_voids_to_add(self):

        root_voids = self.agent.meta['root-mesh-voids']
        alpha_voids = self.cache['alpha-voids-hanging']

        to_extend = np.isin(
            alpha_voids['pivots'], root_voids['pivots'], invert=True
        )

        triangs = alpha_voids['triangs'][to_extend, :]

        return np.copy(
            triangs[:, [1, 0, 2]], order='C'
        )

    def find_alpha_voids_hanging(self):

        inds = self.find_alpha_voids_rank_one()
        mask = self.find_alpha_voids_not_on_edge()

        mask = mask[inds]
        inds = inds[mask]

        triangs = self.triangs_alpha[inds, :]

        voids_meta = {
            'triangs': triangs, 'pivots': triangs[:, 2]
        }

        self.cache['alpha-voids-hanging'] = voids_meta
        return voids_meta

    def find_alpha_voids_rank_one(self):

        data = tables.maptable(
            self.meta['triangs-alpha-codes']
        )

        _, inds, _ = data.atrank(1)
        return inds

    def find_alpha_voids_not_on_edge(self):

        pivots = self.triangs_alpha[:, 2]
        root_edge = self.agent.meta['root-mesh-edge']

        not_on_edge = np.isin(
            pivots, root_edge.nodnums_unique, invert=True
        )

        return not_on_edge


class MakerMeshGamma(SubSupTriAgent):
    """Makes the gamma-mesh from the beta-triangs.

    - Removes super-bodies of the supertriu.
    - Removes root-voids that were doubled.

    - Adds new points.
    - Adds new triangs.

    """

    @property
    def root_points(self):
        return self.agent.suptri.mesh.points

    @property
    def root_mesh_copy(self):
        return self.agent.suptri.mesh.twin()

    @property
    def triangs_alpha(self):
        return self.agent.cache['triangs-alpha']

    @property
    def triangs_beta(self):
        return self.agent.cache['triangs-beta']

    def get_mesh(self):

        _ = self.set_mesher_cache()
        _ = self.run_mesher()

        return _

    def run_mesher(self):

        mesh = self.root_mesh_copy

        mesh = mesh.deltriangs(
            *self.cache['trinums-to-del']
        )

        mesh = mesh.add_points(
            self.cache['points-to-add']
        )

        mesh = mesh.add_triangs(
            self.cache['triangs-to-add']
        )

        return mesh

    def set_mesher_cache(self):

        self.cache |= {
            'trinums-to-del': self.make_trinums_to_del(),
            'triangs-to-add': self.make_triangs_to_add(),
            'points-to-add': self.make_points_to_add()
        }

        return True

    def make_trinums_to_del(self):

        supbodies = self.agent.suptri.supbodies
        garbage_voids = self.triangs_beta['extra-voids']['trinums-to-del']

        data = [
            supbodies.flatten(), garbage_voids
        ]

        return np.hstack(data)

    def make_triangs_to_add(self):

        inner_voids = self.triangs_beta['inner-voids']
        inner_cores = self.triangs_beta['inner-cores']
        inner_sides = self.triangs_beta['inner-sides']
        outer_voids = self.triangs_beta['extra-voids']['triangs-to-add']

        return np.vstack(
            [inner_voids, inner_cores, inner_sides, outer_voids]
        )

    def make_inner_side_triangs(self):

        ker = np.hsplit(
            self.triangs_beta['inner-cores'], 3
        )

        sup = np.hsplit(
            self.agent.suptri.supmesh.triangs, 3
        )

        side1 = _stack_cols(sup[0], ker[0], ker[2])
        side2 = _stack_cols(sup[1], ker[1], ker[0])
        side3 = _stack_cols(sup[2], ker[2], ker[1])

        return np.vstack(
            [side1, side2, side3]
        )

    def make_points_to_add(self):
        return self.root_points[
            self.triangs_alpha[:, 2]
        ]


class MakerMeshRelease(SubSupTriAgent):
    """Makes the mass-mesh from the gamma-mesh.

    - Removes ghost nodes.
    - Inserts metadata.

    """

    @property
    def mesh_gamma(self):
        return self.agent.cache['mesh-gamma']

    @property
    def triangs_alpha(self):
        return self.agent.cache['triangs-alpha']

    @property
    def root_nodes_range(self):
        return np.arange(
            self.agent.suptri.mesh.npoints
        )

    def get_mesh(self):

        mesh = self.del_ghosts()
        mesh = self.add_meta(mesh)

        return mesh

    def del_ghosts(self):

        new_mesh = self.mesh_gamma.delghosts()

        self.cache['nodes-gamma2mass'] = new_mesh.meta['old-nodes-numbers']
        return new_mesh.twin()

    def add_meta(self, mesh):

        root2gamma = np.hstack(
            [self.root_nodes_range, self.triangs_alpha[:, 2]]
        )

        gamma2mass = self.cache['nodes-gamma2mass']

        mesh.meta = {
            'nodes-root2gamma': root2gamma,
            'nodes-gamma2mass': gamma2mass
        }

        return mesh


def _stack_cols(*cols):
    return np.hstack(cols)
