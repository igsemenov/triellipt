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


class MSHReader:
    """Reader of `.msh` meshes.
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

    def with_root(self, path):
        """Assigns the root directory to the reader.

        Parameters
        ----------
        path : str
            Absolute path to the root directory.

        Returns
        -------
        MSHReader
            New reader.

        """
        return self.__class__(path)

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

        mesh_dict = self.read_mesh_dict(file_name)
        trimesh = self.from_mesh_dict(mesh_dict)

        return trimesh

    def from_mesh_dict(self, mesh_dict):
        return trimesh_.TriMesh.from_mesh_dict(mesh_dict)

    def read_mesh_dict(self, file_name) -> dict:

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
