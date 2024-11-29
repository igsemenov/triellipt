# -*- coding: utf-8 -*-
"""Mesh coarsening.
"""
import numpy as np
from scipy import sparse as sp
from triellipt.fem import femunit
from triellipt.amr import supclean
from triellipt.amr import massmesh
from triellipt.amr import utils_


def coarsen_mesh(mesh, trinums):
    _ = get_coarsener(mesh, trinums)
    return _.release_mesh()


def get_target_suptri(mesh, trinums):
    _ = get_coarsener(mesh, trinums)
    return _.make_target_suptri()


def get_mesh_alpha(mesh, trinums):
    _ = get_coarsener(mesh, trinums)
    return _.make_mesh_alpha()


def get_mesh_beta(mesh, trinums):
    _ = get_coarsener(mesh, trinums)
    return _.make_mesh_beta()


def get_mesh_gamma(mesh, trinums):
    _ = get_coarsener(mesh, trinums)
    return _.make_mesh_gamma()


def get_coarsener(mesh, trinums):
    return MeshCoarsener.from_mesh(mesh).with_trinums(trinums)


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
        return {}

    @property
    def root_mesh_copy(self):
        return self.mesh.twin()


class MeshCoarsener(MeshAgent):
    """Coarsens the target triangles.
    """

    def fetch_meta(self):
        return {
            'root-mesh-edge': self.mesh.meshedge()
        }

    def with_trinums(self, trinums):
        _ = self.trinums_filter.get_filtered(trinums)
        self.meta['trinums-to-coarsen'] = _
        return self

    @property
    def trinums_filter(self):
        return FilterTrinums.from_agent(self)

    @property
    def target_trinums(self):
        return self.meta['trinums-to-coarsen']

    def release_mesh(self):
        """Creates the coarsened mesh.
        """

        mesh = self.make_mesh_gamma()

        if mesh is None:
            return self.root_mesh_copy

        mesh = mesh.twin()

        mesh.meta = {
            'data-collect': self.make_data_collector()
        }

        return mesh

    def make_mesh_gamma(self):

        mesh_beta = self.make_mesh_beta()

        if mesh_beta is None:
            return None

        _ = self.maker_mesh_gamma.get_mesh()
        self.cache['mesh-gamma'] = _
        return _

    def make_mesh_beta(self):
        """Removes twin voids.
        """

        mesh_alpha = self.make_mesh_alpha()

        if mesh_alpha is None:
            return None

        _ = self.maker_mesh_beta.get_mesh()
        self.cache['mesh-beta'] = _
        return _

    def make_mesh_alpha(self):

        suptri = self.make_target_suptri()

        if suptri is None:
            return None

        _ = self.maker_mesh_alpha.get_mesh()
        self.cache['mesh-alpha'] = _
        return _

    def make_target_suptri(self):
        _ = self.maker_target_suptri.get_suptri()
        self.cache['target-suptri'] = _
        return _

    def make_data_collector(self):
        return self.maker_data_collect.get_collector()

    @property
    def maker_target_suptri(self):
        return MakerTargetSuptri.from_agent(self)

    @property
    def maker_mesh_alpha(self):
        return MakerMeshAlpha.from_agent(self)

    @property
    def maker_mesh_beta(self):
        return MakerMeshBeta.from_agent(self)

    @property
    def maker_mesh_gamma(self):
        return MakerMeshGamma.from_agent(self)

    @property
    def maker_data_collect(self):
        return MakerDataCollect.from_agent(self)


class MeshSubAgent:
    """Operator on a mesh agent.
    """

    def __init__(self, agent):
        self.agent = agent
        self.meta = self.fetch_meta()
        self.cache = self.push_cache()

    @classmethod
    def from_agent(cls, agent):
        return cls(agent)

    def fetch_meta(self):
        return {}

    def push_cache(self):
        return {}


class MakerTargetSuptri(MeshSubAgent):
    """Makes the target supertriu.
    """

    def fetch_meta(self):
        return {
            'suptri-prime': self.agent.mesh.supertriu()
        }

    @property
    def suptri_prime(self):
        return self.meta['suptri-prime']

    def get_suptri(self):

        if self.agent.target_trinums.size == 0:
            return None

        suptri = self.fetch_target_suptri()
        suptri = self.clean_target_suptri(suptri)

        if suptri.size == 0:
            return None
        return suptri

    def fetch_target_suptri(self):
        return self.suptri_prime.atcores(*self.agent.target_trinums)

    def clean_target_suptri(self, suptri):

        suptri = self.clean_suptri_overlaps(suptri)
        suptri = self.clean_suptri_not_aligned(suptri)

        return suptri

    def clean_suptri_overlaps(self, suptri):
        return supclean.clean_overlaps(suptri)

    def clean_suptri_not_aligned(self, suptri):
        return supclean.clean_not_aligned(suptri)


class MakerMeshAlpha(MeshSubAgent):
    """Makes the alpha-mesh from the root-mesh.

    - Removes super-bodies of the coarse triangles.
    - Inserts new voids on the sides of coarse triangles.

    """

    @property
    def target_suptri(self):
        return self.agent.cache['target-suptri']

    def get_mesh(self):

        self.set_mesher_cache()
        self.run_mesher()

        return self.cache['mesh-alpha']

    def set_mesher_cache(self):

        self.cache |= {
            'triangs-to-add': self.triangs_to_add(),
            'trinums-to-del': self.trinums_to_del()
        }

    def run_mesher(self):

        mesh = self.agent.root_mesh_copy

        trinums_to_del = self.cache['trinums-to-del']
        triangs_to_add = self.cache['triangs-to-add']

        mesh = mesh.deltriangs(*trinums_to_del)
        mesh = mesh.add_triangs(triangs_to_add)

        self.cache['mesh-alpha'] = mesh

    def trinums_to_del(self):
        return self.target_suptri.supbodies.flatten()

    def triangs_to_add(self):
        return self.make_new_voids()

    def make_new_voids(self):

        ker = np.hsplit(
            self.target_suptri.kermesh.triangs, 3
        )

        sup = np.hsplit(
            self.target_suptri.supmesh.triangs, 3
        )

        void_0 = _stack_cols(sup[1], sup[0], ker[1])
        void_1 = _stack_cols(sup[2], sup[1], ker[2])
        void_2 = _stack_cols(sup[0], sup[2], ker[0])

        new_voids = np.vstack(
            [void_0, void_1, void_2]
        )

        self.cache['new-voids'] = new_voids
        return new_voids


class MakerMeshBeta(MeshSubAgent):
    """Makes the beta-mesh from the alpha-mesh.

    - Removes twin voids.
    - Removes voids touching the root-mesh edge.

    """

    def get_mesh(self):

        mesh = self.agent.cache['mesh-alpha']

        mesh = self.clean_twin_voids(mesh)
        mesh = self.clean_voids_on_edge(mesh)

        return mesh

    def clean_twin_voids(self, mesh):
        return utils_.clean_twin_voids(mesh)

    def clean_voids_on_edge(self, mesh):
        return utils_.clean_voids_on_edge(
            mesh, self.agent.meta['root-mesh-edge']
        )


class MakerMeshGamma(MeshSubAgent):
    """Makes the gamma-mesh from the beta-mesh.

    - Adds coarse triangles into the mesh.
    - Removes ghost nodes.

    """

    def fetch_meta(self):
        return {
            'coarse-triangs': self.fetch_coarse_triangs()
        }

    def fetch_coarse_triangs(self):
        return self.agent.cache['target-suptri'].supmesh.triangs

    def push_cache(self):
        return {
            'mesh-pre-gamma': self.mesh_beta.twin()
        }

    @property
    def mesh_beta(self):
        return self.agent.cache['mesh-beta']

    def get_mesh(self):

        mesh = self.cache['mesh-pre-gamma']

        mesh = self.add_coarse_triangs(mesh)
        mesh = self.del_ghost_nodes(mesh)

        return mesh

    def add_coarse_triangs(self, mesh):
        return mesh.add_triangs(
            self.meta['coarse-triangs']
        )

    def del_ghost_nodes(self, mesh):

        new_mesh = mesh.delghosts()

        new_mesh.meta = {
            'nodes-beta2gamma': new_mesh.meta['old-nodes-numbers']
        }

        return new_mesh


class MakerDataCollect(MeshSubAgent):
    """Makes the data-collector.
    """

    @property
    def mesh_gamma(self):
        return self.agent.cache['mesh-gamma']

    def get_collector(self):

        _ = self.make_mass_unit()
        _ = self.make_collector()

        return _

    def make_collector(self):

        inds = self.get_data_inds()
        oprs = self.get_mass_oprs()

        meta = {
            **inds, **oprs
        }

        return DataCollector(
            self.agent.mesh, meta
        )

    def get_data_inds(self):
        return {
            'root2mass': self.get_root2mass(),
            'mass2unit': self.get_mass2unit(),
            'unit2mass': self.get_unit2mass(),
            'root2data': self.get_root2data()
        }

    def get_root2mass(self):

        mass_mesh = self.cache['mass-mesh']

        root2gamma = mass_mesh.meta['nodes-root2gamma']
        gamma2mass = mass_mesh.meta['nodes-gamma2mass']

        return np.copy(
            root2gamma[gamma2mass], order='C'
        )

    def get_mass2unit(self):
        return np.copy(
            self.cache['mass-unit'].perm.perm, order='C'
        )

    def get_unit2mass(self):
        return np.copy(
            self.cache['mass-unit'].perm.perm_inv, order='C'
        )

    def get_root2data(self):
        return np.copy(
            self.mesh_gamma.meta['nodes-beta2gamma'], order='C'
        )

    def get_mass_oprs(self):

        unit = self.cache['mass-unit']

        mass_fem = unit.massopr(
            is_lumped=True, add_constr=False
        )

        mass_amr = unit.massopr(
            is_lumped=True, add_constr=True
        )

        mass_mat = mass_fem.body.copy()
        mass_amr = mass_amr.body.tocsc()

        return {
            'mass-mat': mass_mat,
            'mass-amr': mass_amr
        }

    def make_mass_unit(self):

        mesh = self.make_mass_mesh()

        unit = femunit.getunit(mesh)
        self.cache['mass-unit'] = unit

        return unit

    def make_mass_mesh(self):
        mesh = massmesh.get_massmesh(self.root_suptri)
        self.cache['mass-mesh'] = mesh
        return mesh

    @property
    def root_suptri(self):
        return self.agent.cache['target-suptri']


class DataCollector:
    """Data-collector.
    """

    def __init__(self, mesh, meta):
        self.mesh = mesh
        self.meta = meta
        self.cache = {}

    def __call__(self, data):
        return self.collect(data)

    @property
    def root2mass(self):
        """Images of the source data on the mass-mesh.
        """
        return self.meta['root2mass']

    @property
    def root2data(self):
        """Images of the source data on the coarse mesh.
        """
        return self.meta['root2data']

    @property
    def mass2unit(self):
        return self.meta['mass2unit']

    @property
    def unit2mass(self):
        return self.meta['unit2mass']

    @property
    def massmat(self):
        return self.meta['mass-mat']

    @property
    def massamr(self):
        return self.meta['mass-amr']

    @property
    def massinv(self):
        if 'mass-inv' in self.cache:
            return self.cache['mass-inv']
        self.cache['mass-inv'] = sp.linalg.splu(self.massamr)
        return self.massinv

    def collect(self, data):

        data_fem = data[self.root2mass][self.mass2unit]

        mass_fem = self.massmat @ data_fem
        data_fem = self.massinv.solve(mass_fem)

        data_new = np.zeros_like(
            self.root2mass, dtype=float
        )

        data_new[self.root2mass] = data_fem[self.unit2mass]

        return np.copy(
            data_new[self.root2data], order='C'
        )


class FilterTrinums(MeshSubAgent):
    """Filters the target trinums.
    """

    def get_filtered(self, trinums):

        trinums = np.array(trinums)

        trinums = self.filter_are_unique(trinums)
        trinums = self.filter_are_in_mesh(trinums)

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

    @property
    def mesh_size(self):
        return self.agent.mesh.ntriangs


def _stack_cols(*cols):
    return np.hstack(cols)
