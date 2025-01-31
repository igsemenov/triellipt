# -*- coding: utf-8 -*-
"""Tests the mesh remunerator.
"""
import unittest
from triellipt import mesher


class TestRenum(unittest.TestCase):

    PERM = [
        0, 3, 6, 7, 8, 5, 2, 1, 4
    ]

    @classmethod
    def setUpClass(cls):
        cls.MESH = cls.mesh().renumed(cls.PERM)

    @classmethod
    def mesh(cls):
        return mesher.trigrid(3, 3, 'west-snake')

    def test_renumed(self):
        assert self.loop_nodes.tolist() == list(range(8))

    @property
    def loop_nodes(self):
        return self.mesh_edge.getloops()[0].nodnums1

    @property
    def mesh_edge(self):
        return self.MESH.meshedge()


class TestShuffle(unittest.TestCase):

    def test_edge_primary(self):
        assert self.edge_primary.tolist() == [1, 2, 5, 6]

    def test_edge_aligned(self):
        assert self.edge_aligned.tolist() == [0, 1, 2, 3]

    @property
    def mesh_primary(self):
        return mesher.trigrid(3, 3, 'west-snake')

    @property
    def mesh_aligned(self):
        return self.mesh_primary.shuffled(
            (1, 2, 5, 6, 0, 3, 4, 7)
        )

    @property
    def edge_aligned(self):
        return self.mesh_aligned.meshedge().trinums_unique

    @property
    def edge_primary(self):
        return self.mesh_primary.meshedge().trinums_unique


if __name__ == '__main__':
    unittest.main()
