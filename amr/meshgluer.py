# -*- coding: utf-8 -*
"""Joiner of two meshes.
"""
import numpy as np


def join_meshes(mesh1, mesh2, tol=None):
    """Join the meshes along a shared boundary, if available.

    Parameters
    ----------
    mesh1 : TriMesh
        1-st input mesh.
    mesh2 : TriMesh
        2-nd input mesh.
    tol : int = None
        Optional absolute tolerance for detecting nearby points.

    Returns
    -------
    TriMesh | None
        New mesh or None, if failed.

    """
    return DomainsGluer(mesh1, mesh2, tol).make_new_mesh()


class DomainsGluer:
    """Gluer of two meshes.
    """

    TOL = 9

    def __init__(self, mesh1, mesh2, tol):

        self.mesh1 = mesh1
        self.mesh2 = mesh2
        self.tol = tol or self.TOL

        self.meta = self.fetch_meta()
        self.cache = {}

    def fetch_meta(self):
        return {
            'edge1': self.mesh1.meshedge(),
            'edge2': self.mesh2.meshedge()
        }

    @property
    def edge1(self):
        return self.meta['edge1']

    @property
    def edge2(self):
        return self.meta['edge2']

    @property
    def nodes1(self):
        return self.meta['edge1'].nodes_complex

    @property
    def nodes2(self):
        return self.meta['edge2'].nodes_complex

    @property
    def numshift(self):
        return self.mesh1.npoints

    def make_new_mesh(self):

        state = self.find_twins()

        if state is False:
            return None

        data = self.make_data()
        mesh = self.push_mesh(data)

        return mesh.delghosts()

    def push_mesh(self, data):
        return self.mesh1.from_data(
            data['points'], data['triangs']
        )

    def make_data(self):
        return {
            'points': self.make_new_points(),
            'triangs': self.make_new_triangs()
        }

    def make_new_points(self):
        return np.hstack(
            [self.mesh1.points, self.mesh2.points]
        )

    def make_new_triangs(self):

        recinds = self.make_reconnector()
        triangs = self.make_triangs_long()

        return np.copy(
            recinds[triangs], order='C'
        )

    def make_triangs_long(self):
        return np.vstack(
            [self.mesh1.triangs, self.mesh2.triangs + self.numshift]
        )

    def make_reconnector(self):

        indrange = self.make_range()
        i_1, i_2 = self.take_twins()

        indrange[i_1] = i_2 + self.numshift
        return indrange

    def make_range(self):
        return np.hstack([
            np.arange(self.mesh1.npoints),
            np.arange(self.mesh2.npoints) + self.numshift
        ])

    def take_twins(self):
        return self.cache['twins']

    def find_twins(self):

        nodes1 = self.nodes1.round(self.tol)
        nodes2 = self.nodes2.round(self.tol)

        out, i_1, i_2 = np.intersect1d(
            nodes1, nodes2, return_indices=True
        )

        state = out.size > 0

        twins = [
            self.edge1.nodnums1[i_1],
            self.edge2.nodnums1[i_2]
        ]

        self.cache = {
            'twins': twins,
            'state': state
        }

        return state
