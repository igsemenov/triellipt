# -*- coding: utf-8 -*-
"""Mass collector (solid).

- No voids in meshes.
- Uniform mesh refinement.

"""
import numpy as np
from scipy import sparse as sp
from triellipt.fem import femoprs


def get_collector_structed(mesh1, mesh2):
    _ = MakeCollectStructed(mesh1, mesh2)
    return _.get_collector()


def get_collector_scaled(mesh1, mesh2):
    _ = MakeCollectScaled(mesh1, mesh2)
    return _.get_collector()


class MassCollector:
    """Mass collector object.
    """

    def __init__(self):
        self.mass_diag_parent = None
        self.collector_matrix = None
        self.meta = {}

    @property
    def matrix(self):
        return self.collector_matrix

    @property
    def massdiag(self):
        return self.mass_diag_parent

    def collect_mass(self, slave_data):
        """Collects the mass from the slave mesh.

        Parameters
        ----------
        slave_data : flat-float-array
            Nodal data defined on the slave mesh.

        Returns
        -------
        flat-float-array
            Total mass over the master (parent) nodes.

        """
        return self.matrix @ slave_data

    def get_master_data(self, slave_data):
        """Computes nodal data over the master nodes.

        Parameters
        ----------
        slave_data : flat-float-array
            Nodal data defined on the slave mesh.

        Returns
        -------
        flat-float-array
            Nodal data over the master (parent) nodes.

        """
        nodal_mass = self.collect_mass(slave_data)
        return nodal_mass / self.massdiag


class MeshDuetAgent:
    """Operator on a duet of meshes.

    Attributes
    ----------
    mesh1 : TriMesh
        Parent triangle mesh.
    mesh2 : TriMesh
        Refined triangle mesh.

    """

    def __init__(self, mesh1, mesh2):

        self.mesh1 = mesh1
        self.mesh2 = mesh2

        self.meta = self.fetch_meta()
        self.cache = {}

    @classmethod
    def from_mesh_duet(cls, mesh1, mesh2):
        return cls(mesh1, mesh2)

    def fetch_meta(self):
        return {
            'mesh2-suptriu': self.mesh2.supertriu()
        }


class CollectorMaker(MeshDuetAgent):
    """Maker of a mass-collector.
    """

    def get_collector(self):

        obj = MassCollector()

        obj.collector_matrix = self.get_collector_matrix()
        obj.mass_diag_parent = self.get_mass_diag_parent()

        obj.meta = self.get_meta()
        return obj

    def get_collector_matrix(self):

        _ = self.dump_duet_image()
        _ = self.push_matrix()

        return _

    def push_matrix(self):

        v_stream = self.make_v_stream()
        ij_stream = self.make_ij_stream()

        as_coo = sp.coo_array(
            (v_stream, ij_stream), shape=self.collector_shape
        )

        as_csr = as_coo.tocsr()
        return as_csr

    def dump_duet_image(self):
        _ = self.duet_imager.get_duet_image()
        self.cache['duet-image'] = _
        return _

    def make_ij_stream(self):

        hosts = np.hsplit(self.hosts_mesh.triangs, 3)
        cores = np.hsplit(self.cores_mesh.triangs, 3)

        i_stream_0 = _stack_cols(hosts[0], hosts[0], hosts[0])
        i_stream_1 = _stack_cols(hosts[1], hosts[1], hosts[1])
        i_stream_2 = _stack_cols(hosts[2], hosts[2], hosts[2])

        j_stream_0 = _stack_cols(hosts[0], cores[0], cores[1])
        j_stream_1 = _stack_cols(hosts[1], cores[1], cores[2])
        j_stream_2 = _stack_cols(hosts[2], cores[2], cores[0])

        i_stream = _stack_cols(
            i_stream_0, i_stream_1, i_stream_2
        )

        j_stream = _stack_cols(
            j_stream_0, j_stream_1, j_stream_2
        )

        i_stream = i_stream.flatten()
        j_stream = j_stream.flatten()

        return (i_stream, j_stream)

    def make_v_stream(self):

        mass_coeffs_proxy = np.ones(
            (self.hosts_mesh.ntriangs, 9)
        )

        return mass_coeffs_proxy.flatten()

    def get_mass_diag_parent(self):
        return self.mesh1.nodesmap().nodes_ranks

    def get_meta(self):
        """Returns collector metadata.
        """
        return self.get_mesh_meta() | self.get_collector_mode()

    def get_mesh_meta(self):
        return {
            'mesh-1': _mesh_meta(self.mesh1),
            'mesh-2': _mesh_meta(self.mesh2)
        }

    def get_collector_mode(self):
        return {
            'mode': ''
        }

    @property
    def duet_imager(self):
        return DuetImager.from_agent(self)

    @property
    def hosts_mesh(self):
        return self.cache['duet-image']['hosts']

    @property
    def cores_mesh(self):
        return self.cache['duet-image']['cores']

    @property
    def collector_shape(self):
        return (
            self.mesh1.npoints, self.mesh2.npoints
        )


class MakeCollectStructed(CollectorMaker):
    """Makes a mass collector for meshes with equal triangles.
    """

    def make_v_stream(self):

        mass_stream = np.tile(
            self.mass_node_coeffs, (self.hosts_mesh.ntriangs, 3)
        )

        return mass_stream.flatten()

    def get_mass_diag_parent(self):
        return self.mesh1.nodesmap().nodes_ranks / 3.

    @property
    def mass_node_coeffs(self):
        return np.r_[2., 3., 3.] / 24.

    def get_collector_mode(self):
        return {
            'mode': 'structed'
        }


class MakeCollectScaled(CollectorMaker):
    """Makes a mass collector for general meshes.
    """

    def make_v_stream(self):

        mass_stream = np.tile(
            self.mass_node_coeffs, (self.hosts_mesh.ntriangs, 3)
        )

        host_areas = self.get_host_mesh_areas()

        mass_stream = mass_stream * host_areas[..., None]
        return mass_stream.flatten()

    def get_mass_diag_parent(self):

        areas = self.get_parent_mesh_areas()
        nodes = self.get_parent_nodesmap()

        nodes_areas = np.add.reduceat(
            areas[nodes.trinums], nodes.data.bins_reduce
        )

        return nodes_areas / 3.

    def get_host_mesh_areas(self):
        return _mesh_areas(self.hosts_mesh)

    def get_parent_mesh_areas(self):
        return _mesh_areas(self.mesh1)

    def get_parent_nodesmap(self):
        return self.mesh1.nodesmap()

    @property
    def mass_node_coeffs(self):
        return np.r_[2., 3., 3.] / 24.

    def get_collector_mode(self):
        return {
            'mode': 'unstructed'
        }


class DuetSubAgent:
    """Operator on a mesh-duet agent.
    """

    def __init__(self, agent):
        self.agent = agent
        self.cache = {}

    @classmethod
    def from_agent(cls, agent):
        return cls(agent)


class DuetImager(DuetSubAgent):
    """Makes a duet image as a hosts-cores pair of meshes. 
    """

    @classmethod
    def from_mesh_duet(cls, mesh1, mesh2):
        return cls(
            MeshDuetAgent.from_mesh_duet(mesh1, mesh2)
        )

    def get_duet_image(self):

        _ = self.make_suptriu2_host_triangs_only()
        _ = self.push_duet_image()

        return _

    def push_duet_image(self):
        return {
            'hosts': self.suptriu2_hosts_only.supmesh,
            'cores': self.suptriu2_hosts_only.kermesh
        }

    def make_suptriu2_host_triangs_only(self):

        nodesmap = self.suptriu2.supmesh.nodesmap()

        dataspan_host_nodes = np.sum(
            nodesmap.nodes_ranks[self.range_hostnodes]
        )

        host_triangs_in_suptri2 = np.unique(
            nodesmap.trinums[:dataspan_host_nodes]
        )

        _ = self.suptriu2.subtriu(*host_triangs_in_suptri2)
        self.cache['mesh2-suptriu-hosts-only'] = _
        return _

    @property
    def suptriu2(self):
        """Supertriu for mesh2.
        """
        return self.agent.meta['mesh2-suptriu']

    @property
    def suptriu2_hosts_only(self):
        """Supertriu for mesh2 with host triangles only.
        """
        return self.cache['mesh2-suptriu-hosts-only']

    @property
    def range_hostnodes(self):
        """Range of host (top) nodes in all numberings.
        """
        return np.s_[:self.agent.mesh1.npoints]


def _stack_cols(*cols):
    return np.hstack(cols)


def _mesh_meta(mesh):
    return {
        'npoints': mesh.npoints,
        'ntriangs': mesh.ntriangs
    }


def _mesh_areas(mesh):
    return 0.5 * femoprs.mesh_metric(mesh)['jacobis']
