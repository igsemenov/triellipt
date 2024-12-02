# -*- coding: utf-8 -*-
"""Making a mesh skeleton.
"""
from abc import abstractmethod, ABC
import numpy as np
from triellipt.fem import supvoids
from triellipt.utils import tables
from triellipt.utils import pairs


def getskeleton(mesh):
    """Converts a mesh into a mesh skeleton.
    """
    return MeshSkeleton.from_mesh(mesh)


class MeshSkeleton:
    """Mesh in skeleton format.
    """

    def __init__(self, body, joints):
        self.body = body
        self.joints = joints

    @property
    def wests(self):
        return self.joints.get('wests')

    @property
    def easts(self):
        return self.joints.get('easts')

    @property
    def cores(self):
        return self.joints.get('cores')

    @property
    def voids(self):
        return self.joints.get('voids')

    @property
    def hasjoints(self):
        return bool(self.joints)

    @classmethod
    def from_mesh(cls, mesh):

        if mesh.hasvoids():
            return MeshConverter(mesh).get_skeleton()

        return cls(
            mesh, {}
        )


class MeshConverter:
    """Converts a mesh into a mesh skeleton.
    """

    def __init__(self, mesh):
        self.mesh = mesh

    def get_skeleton(self):

        joints_data = self.make_joints_data()
        mesh_skeleton = self.from_joints_data(joints_data)

        return mesh_skeleton

    def make_joints_data(self):
        return self.joints_collector.get_joints_data()

    def from_joints_data(self, joints_data):

        body = self.make_body(joints_data)
        joints = self.make_joints(joints_data)

        return MeshSkeleton(body, joints)

    def make_body(self, joints_data):

        trinums_to_remove = [
            o['trinums'] for o in joints_data.values()
        ]

        return self.mesh.deltriangs(
            *np.hstack(trinums_to_remove)
        )

    def make_joints(self, joints_data):

        triangs = {
            k: v['triangs'] for k, v in joints_data.items()
        }

        return {
            k: self.mesh.update_triangs(v) for k, v in triangs.items()
        }

    @property
    def joints_collector(self):
        return JointsCollector.from_mesh(self.mesh)


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh):

        self.mesh = mesh
        self.suptri_prime = None
        self.suptri_voids = None

        self.add_suptri_items()

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    def add_suptri_items(self):

        self.add_suptri_prime()
        self.add_suptri_voids()

    def add_suptri_prime(self):
        self.suptri_prime = self.mesh_supertriu

    def add_suptri_voids(self):
        self.add_suptri_voids_proxy()

    def add_suptri_voids_prime(self):
        self.suptri_voids = self.suptri_prime.supvoids()

    def add_suptri_voids_proxy(self):
        self.suptri_voids = supvoids.get_sup_voids(self.mesh)

    @property
    def mesh_supertriu(self):
        return self.mesh.supertriu()

    @property
    def pivots(self):
        return self.mesh.triangs[
            self.suptri_voids.trinums, 2
        ]


class JointsCollector(MeshAgent):
    """Collects joints from the mesh.
    """

    def get_joints_data(self):
        return {
            'wests': self.get_wests(),
            'easts': self.get_easts(),
            'cores': self.get_cores(),
            'voids': self.get_voids()
        }

    def get_wests(self):
        return WestsFinder(self).get_data()

    def get_easts(self):
        return EastsFinder(self).get_data()

    def get_cores(self):
        return CoresFinder(self).get_data()

    def get_voids(self):
        return {
            'trinums': self.suptri_voids.trinums,
            'triangs': self.suptri_voids.kermesh.triangs
        }


class SubAgent(ABC):

    def __init__(self, agent):
        self.agent = agent

    @classmethod
    def from_agent(cls, agent):
        return cls(agent)

    @property
    def pivots(self):
        return self.agent.pivots


class JointsFinder(SubAgent):
    """Base finder of joints items.
    """

    def get_data(self):

        trinums = self.find_trinums()
        triangs = self.make_triangs(trinums)

        return {
            'trinums': trinums,
            'triangs': triangs
        }

    def make_triangs(self, trinums):

        prime_triangs = self.take_prime_triangs(trinums)
        syncedtriangs = self.sync_prime_triangs(prime_triangs)

        return syncedtriangs

    def take_prime_triangs(self, trinums):
        return self.agent.mesh.triangs[trinums, :]

    def sync_prime_triangs(self, prime_triangs):

        pivslocnums = self.find_pivots_in_triangs(prime_triangs)
        new_triangs = self.roll_pivots_in_triangs(prime_triangs, pivslocnums)

        return new_triangs

    def find_pivots_in_triangs(self, triangs_with_pivots):
        """Returns pivots local numbers in target triangles.
        """

        _, pivs_locnums = np.where(
            triangs_with_pivots == self.pivots[..., None]
        )

        return pivs_locnums

    def roll_pivots_in_triangs(self, triangs_with_pivots, pivots_locnums):

        roller = _make_triplets_table(pivots_locnums)

        return tables.table_image(
            triangs_with_pivots, roller
        )

    @abstractmethod
    def find_trinums(self):
        """Returns triangle numbers for a joint item.
        """


class WestsFinder(JointsFinder):
    """Finds west triangles in joints.
    """

    def find_trinums(self):
        return self.agent.suptri_voids.trinums2


class EastsFinder(JointsFinder):
    """Finds east triangles in joints.
    """

    def find_trinums(self):
        return self.agent.suptri_voids.trinums3


class CoresFinder(JointsFinder):
    """Finds core triangles in joints.
    """

    def find_trinums(self):

        trinums_prime = self.find_trinums_prime()
        trinums_order = self.sort_trinums_prime(trinums_prime)

        return trinums_order

    def find_trinums_prime(self):

        neighbours = self.pair_neighbours()
        voids_ears = self.pair_voids_ears()

        mask_table = np.isin(
            neighbours, voids_ears
        )

        mask_conv = np.any(
            mask_table, axis=1
        )

        trinums_with_voids = self.agent.suptri_prime.trinums[mask_conv]

        return np.setdiff1d(
            trinums_with_voids, self.agent.suptri_voids.trinums
        )

    def sort_trinums_prime(self, trinums_prime):

        if trinums_prime.size == 0:
            return trinums_prime

        triangs_prime = self.agent.mesh.triangs[trinums_prime, :]

        permuter = _sync_table_to_pivots(
            triangs_prime, self.agent.pivots
        )

        return trinums_prime[permuter]

    def pair_voids_ears(self):

        supv = self.agent.suptri_voids

        return pairs.sympaired(
            supv.trinums2, supv.trinums3
        )

    def pair_neighbours(self):

        supt = self.agent.suptri_prime

        trio = _pack_cols(
            supt.trinums1, supt.trinums2, supt.trinums3
        )

        return pairs.paircols(
            trio[:, [0, 1, 2, 0]]
        )


def _sync_table_to_pivots(table, pivots):

    pivots_in_table = table[
        *np.where(np.isin(table, pivots))
    ]

    _, ind1, ind2 = np.intersect1d(
        pivots_in_table, pivots, return_indices=True
    )

    rows_permuter = ind1[np.argsort(ind2)]
    return rows_permuter


def _make_triplets_table(start_indices):

    _ = start_indices

    return _pack_cols(
        _, (_ + 1) % 3, (_ + 2) % 3
    )


def _pack_cols(*cols):
    return np.vstack(cols).T.copy('C')
