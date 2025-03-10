# -*- coding: utf-8 -*-
"""Reads meshes from MSH files.
"""
import os
from triellipt.mshread import mshparser
from triellipt.mshread import mshblocks
from triellipt.trimesh import trimesh_

MSHReaderError = type(
    'MSHReaderError', (Exception,), {}
)


def getreader(path):
    """Creates a reader of Gmsh meshes.

    Parameters
    ----------
    path : str
        Path to the folder with Gmsh meshes.

    Returns
    -------
    MSHReader
        Reader of Gmsh meshes.

    """
    return MSHReader.from_path(path)


class MSHReader:
    """Reader of Gmsh meshes.
    """

    def __init__(self, root_path=None):

        self.root_path = root_path

        if root_path is None:
            return

        try:
            path_is_dir = os.path.isdir(root_path)
        except TypeError as exc:
            raise MSHReaderError(str(exc)) from exc

        if path_is_dir:
            return

        raise MSHReaderError(
            "root directory not found, check the path"
        )

    @classmethod
    def from_path(cls, path):
        return cls(path)

    def listmeshes(self):
        """Returns list of `.msh` files in the root directory.
        """

        if not self.root_path:
            return []

        files = os.listdir(self.root_path)

        return [
            f for f in files if f.endswith('.msh')
        ]

    def read_mesh(self, file_name):
        """Reads a mesh from an `.msh` file.

        Parameters
        ----------
        file_name : str
            Name of the `.msh` file.

        Returns
        -------
        TriMesh
            Mesh object.

        """

        mshdata = self.read_mesh_data(file_name)
        trimesh = self.from_mesh_data(mshdata)

        trimesh = trimesh.delghosts()
        return trimesh

    def from_mesh_data(self, mesh_dict):
        return trimesh_.TriMesh.from_mesh_dict(mesh_dict)

    def read_mesh_data(self, file_name) -> dict:

        try:
            filepath = os.path.join(
                self.root_path, file_name
            )
        except TypeError as exc:
            raise MSHReaderError(str(exc)) from exc

        if not os.path.exists(filepath):
            raise MSHReaderError(
                "file not found, check the path"
            )

        sections = self.parse_into_sections(filepath)
        meshdict = self.fetch_from_sections(sections)

        return meshdict

    def parse_into_sections(self, filepath):
        content = self.read_file(filepath)
        return mshparser.MSHParser.with_content(content).getsections()

    def fetch_from_sections(self, sections: dict):

        nodes = mshblocks.MSHNodes.from_sections(sections)
        elements = mshblocks.MSHElements.from_sections(sections)

        return {
            'nodes': nodes,
            'elements': elements
        }

    def read_file(self, abspath):
        with open(abspath, encoding='utf-8', mode='r') as file:
            return file.read()
