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


if __name__ == '__main__':
    unittest.main()
