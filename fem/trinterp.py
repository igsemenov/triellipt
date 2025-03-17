# -*- coding: utf-8 -*-
"""Interpolation on a mesh.
"""
import numpy as np


def getinterp(mesh, xseeds, yseeds):
    """Creates an interpolator for a mesh.

    Parameters
    ----------
    mesh : TriMesh
        Parent mesh.
    xseeds : float-flat-array
        X-coordinates of interpolation nodes.
    yseeds : float-flat-array
        Y-coordinates of interpolation nodes.

    Returns
    -------
    TriInterp
        Interpolator on a mesh of the FEM unit.

    """

    maker = get_interp_maker(mesh)

    return maker.get_triinterp(
        *align_nodes(xseeds, yseeds)
    )


def align_nodes(xnodes, ynodes):

    size = min(
        xnodes.size, ynodes.size
    )

    xnodes = xnodes[:size].copy()
    ynodes = ynodes[:size].copy()

    return xnodes, ynodes


def get_interp_maker(mesh):
    """Interpolator maker for a mesh.
    """

    if not mesh.hasvoids():
        return TriInterpMaker.from_mesh(mesh)

    mesh_proxy = mesh.delvoids()
    return TriInterpMaker.from_mesh(mesh_proxy)


class TriInterp:
    """Interpolator on a mesh.
    """

    def __init__(self, mesh=None, meta=None):
        self.mesh = mesh
        self.meta = meta

    @property
    def xnodes(self):
        return self.points[0, :]

    @property
    def ynodes(self):
        return self.points[1, :]

    @property
    def points(self):
        return self.meta['points']

    @property
    def triangs(self):
        return self.mesh.triangs[
            self.meta['trinums'], :
        ]

    def __call__(self, nodes_data):
        return self.interp(nodes_data)

    def interp(self, nodes_data):
        return np.sum(
            self.meta['coeffs'] * nodes_data[self.triangs], axis=1
        )

    def feval(self, func):
        return func(self.xnodes, self.ynodes)


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh):
        self.mesh = mesh

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    @property
    def vertices(self):
        return self.mesh.points[self.mesh.triangs]


class TriInterpMaker(MeshAgent):
    """Interpolator maker.
    """

    def get_triinterp(self, xseeds, yseeds):

        find_meta = self.find_points(xseeds, yseeds)
        interpobj = self.from_find_meta(find_meta)

        return interpobj

    def find_points(self, xseeds, yseeds):
        return self.tri_finder.find_points(xseeds, yseeds)

    def from_find_meta(self, find_meta):

        coeffs = self.make_coeffs(find_meta)

        meta = {
            **find_meta, 'coeffs': coeffs
        }

        return TriInterp(self.mesh, meta)

    def make_coeffs(self, find_meta):

        verts = self.take_verts(find_meta)
        nodes = self.take_nodes(find_meta)

        return InterpCoeffs(nodes, verts).get_coeffs()

    def take_nodes(self, find_meta):
        return self.pack_complex(*find_meta['points'])

    def take_verts(self, find_meta):

        verts = self.mesh.points[
            self.mesh.triangs[find_meta['trinums'], :]
        ]

        return np.hsplit(verts, 3)

    def pack_complex(self, xpos, ypos):
        return np.atleast_2d(xpos + 1j * ypos).T.copy('C')

    @property
    def tri_finder(self):
        return TriFinder.from_mesh(self.mesh)


class TriFinder(MeshAgent):
    """Finder of host triangles.
    """

    def __init__(self, mesh):
        super().__init__(mesh)
        self.points = None

    def find_points(self, xseeds, yseeds):
        self.points = np.atleast_2d(xseeds + 1j * yseeds)

        mesh_meta = self.make_mesh_meta()
        find_meta = self.from_mesh_meta(mesh_meta)

        return find_meta

    def make_mesh_meta(self):

        zverts = self.vertices

        zdiffs = np.diff(
            zverts[:, [0, 1, 2, 0]], axis=1
        )

        scales = np.abs(zdiffs) / zdiffs

        return {
            'zverts': zverts,
            'scales': scales
        }

    def from_mesh_meta(self, mesh_meta):

        mask_data = self.make_mask_data(mesh_meta)
        find_meta = self.from_mask_data(mask_data)

        return find_meta

    def make_mask_data(self, mesh_meta):

        meta_long = self.make_mesh_meta_long(mesh_meta)
        mask_data = self.mask_from_meta_long(meta_long)

        return mask_data

    def make_mesh_meta_long(self, mesh_meta):

        zverts = mesh_meta['zverts']
        scales = mesh_meta['scales']

        zverts_long = [
            np.tile(v, self.points.shape) for v in np.hsplit(zverts, 3)
        ]

        scales_long = [
            np.tile(v, self.points.shape) for v in np.hsplit(scales, 3)
        ]

        return {
            'zverts': zverts_long,
            'scales': scales_long
        }

    def mask_from_meta_long(self, meta_long):

        zverts = meta_long['zverts']
        scales = meta_long['scales']

        images = [
            (self.points - z) * s for z, s in zip(zverts, scales)
        ]

        is_ins = [
            (v.imag >= 0.0).astype(int) for v in images
        ]

        mask_data = sum(is_ins).T
        return mask_data

    def from_mask_data(self, mask_data):

        ipos, jpos = np.where(mask_data == 3)

        xy_data = self.make_actual_nodes(ipos)

        return {
            'points': xy_data, 'trinums': jpos
        }

    def make_actual_nodes(self, actual_inds):

        points = self.points.flat[actual_inds]

        return np.vstack(
            [points.real, points.imag]
        )


class InterpCoeffs:
    """Computes interpolation coefficients.

    Parameters
    ----------
    nodes : complex-column
        Interpolation nodes.
    verts : list-of-complex-columns
        Vertices of host triangles.

    """

    def __init__(self, nodes, verts):
        self.nodes = nodes
        self.verts = verts

    def get_coeffs(self):

        tri_areas = self.get_areas_triangs()
        nod_areas = self.get_areas_nodes()

        nod_coefs = [
            area / tri_areas for area in nod_areas
        ]

        return np.vstack(nod_coefs).T.copy('C')

    def get_areas_triangs(self):
        return TriangArea(self.verts_table).get_area()

    def get_areas_nodes(self):
        return list(
            self.gen_areas_nodes()
        )

    def gen_areas_nodes(self):
        yield TriangArea(self.node_triangs_one).get_area()
        yield TriangArea(self.node_triangs_two).get_area()
        yield TriangArea(self.node_triangs_tri).get_area()

    @property
    def node_triangs_one(self):
        return np.hstack(
            [self.nodes, self.verts[1], self.verts[2]]
        )

    @property
    def node_triangs_two(self):
        return np.hstack(
            [self.nodes, self.verts[2], self.verts[0]]
        )

    @property
    def node_triangs_tri(self):
        return np.hstack(
            [self.nodes, self.verts[0], self.verts[1]]
        )

    @property
    def verts_table(self):
        return np.hstack(self.verts)


class TriangArea:
    """Triangles area from complex vertices.
    """

    def __init__(self, verts):

        data = verts
        conj = verts.conj()

        self.data = np.hsplit(data, 3)
        self.conj = np.hsplit(conj, 3)

    def get_area(self):

        minor1 = self.get_minor_one()
        minor2 = self.get_minor_two()
        minor3 = self.get_minor_tri()

        area = np.real(
            0.25 * 1j * (minor1 + minor2 + minor3)
        )

        return area.flatten()

    def get_minor_one(self):
        return self.data[1] * self.conj[2] - self.data[2] * self.conj[1]

    def get_minor_two(self):
        return self.data[2] * self.conj[0] - self.data[0] * self.conj[2]

    def get_minor_tri(self):
        return self.data[0] * self.conj[1] - self.data[1] * self.conj[0]
