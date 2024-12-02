# -*- coding: utf-8 -*-
"""Tests removing mesh mouths.

TEST-MESH:

      3---7---4
     / \ / \ / \
    /   5---6   \
   /     \ /     \
  0-------1-------2

"""
import unittest
import numpy as np
from triellipt import trimesh
from triellipt.trimesh import delmouths_


class TestMesh:
    """Defines test meshes.
    """

    POINTS = np.r_[
        0, 1, 2, 0.5 + 1j, 1.5 + 1j, 0.75 + 0.5j, 1.25 + 0.5j, 1 + 1j
    ]

    TRIANGS_BASES = np.array(
        [[0, 1, 3], [1, 2, 4], [1, 4, 3]]
    )

    TRIANGS_CORES = np.array(
        [[1, 6, 5], [6, 4, 7], [5, 7, 3], [5, 6, 7]]
    )

    TRIANGS_VOIDS = np.array([[1, 4, 6], [1, 3, 5]])

    @property
    def mesh_full(self):
        return self.mesh_core.add_triangs(self.TRIANGS_VOIDS)

    @property
    def mesh_core(self):
        return self.mesh_base.add_triangs(self.TRIANGS_CORES)

    @property
    def mesh_base(self):
        return trimesh.TriMesh.from_data(
            self.POINTS, self.TRIANGS_BASES
        )


class TestAgent(TestMesh):

    MESH_AGENT_NAME = None

    @property
    def agent_base(self):
        return self.mesh_agent(self.mesh_base).enriched()

    @property
    def agent_full(self):
        return self.mesh_agent(self.mesh_full).enriched()

    def mesh_agent(self, mesh):
        return getattr(delmouths_, self.MESH_AGENT_NAME)(mesh)


class TestMeshAgent(TestAgent, unittest.TestCase):

    MESH_AGENT_NAME = 'MeshAgent'

    def test_voids_trinums(self):
        assert self.agent_base.voids_trinums.tolist() == []
        assert self.agent_full.voids_trinums.tolist() == [7, 8]

    def test_voids_triangs(self):
        assert self.agent_base.voids_triangs.tolist() == []
        assert self.agent_full.voids_triangs.tolist() == [
            [1, 4, 6], [1, 3, 5]
        ]


class TestGetCorners(TestAgent, unittest.TestCase):

    MESH_AGENT_NAME = 'GetCorners'

    def test_find_corners(self):
        assert self.agent_base.find_corners().tolist() == []
        assert self.agent_full.find_corners().tolist() == [3]


class TestGetTongues(TestAgent, unittest.TestCase):

    MESH_AGENT_NAME = 'GetTongues'

    def test_find_tongues_no_voids(self):
        assert self.agent_base.find_tongues() is None

    def test_find_tongues_with_voids(self):
        assert self.agent_full.find_tongues().trinums.tolist() == [6]
        assert self.agent_full.find_tongues().locnums_apexes.tolist() == [2]
        assert self.agent_full.find_tongues().nodnums_apexes.tolist() == [7]


class TestCleaner(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.CLEANER = delmouths_.make_cleaner(cls.mesh())

    @classmethod
    def mesh(cls):
        return TestMesh().mesh_full

    def test_pick_voids_to_del(self):
        assert self.CLEANER.pick_voids_to_del().tolist() == [7, 8]

    def test_pick_cores_to_del(self):
        assert self.CLEANER.pick_cores_to_del().tolist() == [6, 3, 4, 5]

    def test_make_new_voids(self):
        assert self.CLEANER.make_new_voids().tolist() == [[3, 4, 7]]

    def test_pick_new_cores(self):
        assert self.CLEANER.pick_new_cores().tolist() == [[1, 4, 3]]

    def test_new_mesh_has_no_ghosts(self):
        assert self.mesh_cleaned.hasghosts() is False

    @property
    def mesh_cleaned(self):
        return self.CLEANER.mesh_cleaned()


if __name__ == '__main__':
    unittest.main()
