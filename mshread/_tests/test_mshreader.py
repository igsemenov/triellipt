# -*- coding: utf-8 -*-
"""Test MSH reader.
"""
import os
import unittest
from triellipt.mshread import MSHReader


def dirpath():
    return os.path.join(
        os.path.dirname(__file__), 'msh'
    )


class TestReader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.READER = MSHReader.from_path(dirpath())

    def test_listmeshes(self):
        assert 'circ-gmsh.msh' in self.READER.listmeshes()
        assert 'rect-gmsh.msh' in self.READER.listmeshes()

    def test_mesh_dict_nodes(self):
        assert self.mesh_circ['nodes'].shape == (42, 2)
        assert self.mesh_circ['nodes'].dtype.name == 'float64'

    def test_mesh_dict_elements(self):
        assert self.mesh_circ['elements'].shape == (64, 3)
        assert self.mesh_circ['elements'].dtype.name == 'int32'

    @property
    def mesh_circ(self):
        return self.READER.read_mesh_data('circ-gmsh.msh')

    @property
    def mesh_rect(self):
        return self.READER.read_mesh_data('rect-gmsh.msh')


if __name__ == '__main__':
    unittest.main()
