# -*- coding: utf-8 -*-
"""Triangle lattice.
"""
import numpy as np
from triellipt.utils import pairs, tables
from triellipt.trimesh import trimesh_


def get_lattice(xsize, ysize, close):

    _ = get_latticer(xsize, ysize)
    mesh = _.get_lattice()

    if close is False:
        return mesh
    return close_lattice(mesh)


def get_latticer(xsize, ysize):
    _ = TriLatticer()
    return _.with_counts(xsize, ysize)


def close_lattice(mesh):
    _ = LatticeCloser.from_mesh(mesh)
    return _.get_lattice_closed()


class _TriLatticer:
    """Root of a lattice maker.
    """

    def __init__(self):
        self.icount = None
        self.jcount = None
        self.cache = {}

    def with_counts(self, icount, jcount):

        self.icount = icount
        self.jcount = jcount

        self.set_nodes_grids()
        return self

    def set_nodes_grids(self):
        self.set_nodes_even_even()
        self.set_nodes_odd_odd()

    def set_nodes_even_even(self):

        ieven, jeven = np.meshgrid(
            np.arange(0, self.icount, 2),
            np.arange(0, self.jcount, 2)
        )

        nodes = pairs.szupaired(ieven, jeven)

        self.cache['nodes-even-even'] = {
            'ieven': ieven, 'jeven': jeven, 'nodes': nodes
        }

    def set_nodes_odd_odd(self):

        iodd, jodd = np.meshgrid(
            np.arange(1, self.icount, 2),
            np.arange(1, self.jcount, 2)
        )

        nodes = pairs.szupaired(iodd, jodd)

        self.cache['nodes-odd-odd'] = {
            'iodd': iodd, 'jodd': jodd, 'nodes': nodes
        }

    def sort_mesh(self, mesh):

        mesh = mesh.renumed(
            self.argsort_points(mesh.points_complex)
        )

        mesh = mesh.shuffled(
            self.argsort_centrs(mesh.centrs_complex)
        )

        return mesh.twin()

    def argsort_points(self, points):
        return self.data_sorter(points)

    def argsort_centrs(self, centrs):
        return self.data_sorter(centrs)

    @property
    def data_sorter(self):
        def argsort_complex(data):
            return np.argsort(data)
        return argsort_complex

    @property
    def nodes_oo(self):
        """Odd-odd paired nodes.
        """
        return self.cache['nodes-odd-odd']['nodes']

    @property
    def nodes_ee(self):
        """Even-even paired nodes.
        """
        return self.cache['nodes-even-even']['nodes']

    @property
    def mesh_count(self):
        return self.nodes_ee.size + self.nodes_oo.size


class TriLatticer(_TriLatticer):
    """Maker of a triangle lattice.
    """

    def get_lattice(self):

        mesh = self.make_mesh()
        mesh = self.sort_mesh(mesh)

        meta = {
            'lattice-size': (self.icount, self.jcount)
        }

        return mesh.add_meta(meta)

    def make_mesh(self):

        # First triangles, then points.
        # Points renumbering depends on triangles.
        triangs = self.make_triangs()
        points = self.make_points()

        return trimesh_.TriMesh.from_data(points, triangs)

    def make_points(self):

        points = self.make_points_array()
        sorter = self.make_points_sorter()

        return points[sorter]

    def make_points_array(self):

        epoints = self.make_points_even_even()
        opoints = self.make_points_odd_odd()

        return np.r_[
            epoints.flat, opoints.flat
        ]

    def make_points_sorter(self):

        nodes_long = np.r_[
            self.nodes_ee.flat, self.nodes_oo.flat
        ]

        if self.cache['triangs-complete']:
            return np.argsort(nodes_long)

        in_mesh = np.isin(
            nodes_long, self.cache['triangs-encoded']
        )

        not_in_mesh = np.logical_not(in_mesh)

        nodes_long[not_in_mesh] = nodes_long[in_mesh].max() + 1
        return np.argsort(nodes_long)

    def make_points_even_even(self):

        ieven = self.cache['nodes-even-even']['ieven']
        jeven = self.cache['nodes-even-even']['jeven']

        return ieven + 1j * jeven * np.sqrt(3)

    def make_points_odd_odd(self):

        iodd = self.cache['nodes-odd-odd']['iodd']
        jodd = self.cache['nodes-odd-odd']['jodd']

        return iodd + 1j * jodd * np.sqrt(3)

    def make_triangs(self):

        triangs_encoded = self.make_triangs_encoded()
        triangs_natural = self.make_triangs_natural(triangs_encoded)

        _ = np.amax(triangs_natural) == self.mesh_count - 1
        self.cache['triangs-complete'] = _

        return triangs_natural

    def make_triangs_encoded(self):

        triangs_eeo = self.make_triangs_eeo()
        triangs_ooe = self.make_triangs_ooe()
        triangs_eoo = self.make_triangs_eoo()
        triangs_oee = self.make_triangs_oee()

        triangs = np.vstack(
            [triangs_eeo, triangs_eoo, triangs_ooe, triangs_oee]
        )

        triangs = triangs[:, [2, 1, 0]]
        triangs = triangs.astype(int)

        self.cache['triangs-encoded'] = triangs
        return triangs

    def make_triangs_natural(self, triangs_encoded):
        return tables.norm_table(triangs_encoded)

    def make_triangs_eeo(self):

        onodes = self.nodes_oo
        enodes = self.nodes_ee

        isize = enodes.shape[1] - 1
        jsize = onodes.shape[0]

        verts0 = enodes[:jsize, 1::]
        verts1 = enodes[:jsize, :-1]
        verts2 = onodes[:jsize, :isize]

        return _make_triangs(verts0, verts1, verts2)

    def make_triangs_ooe(self):

        onodes = self.nodes_oo
        enodes = self.nodes_ee

        isize = onodes.shape[1]
        jsize = enodes.shape[0]

        verts0 = onodes[:jsize-1, 1::]
        verts1 = onodes[:jsize-1, :-1]
        verts2 = enodes[1:jsize, 1:isize]

        return _make_triangs(verts0, verts1, verts2)

    def make_triangs_eoo(self):

        onodes = self.nodes_oo
        enodes = self.nodes_ee

        isize = onodes.shape[1]
        jsize = onodes.shape[0]

        verts0 = enodes[:jsize, 1:isize]

        verts1 = onodes[:jsize, :-1]
        verts2 = onodes[:jsize, 1::]

        return _make_triangs(verts0, verts1, verts2)

    def make_triangs_oee(self):

        onodes = self.nodes_oo
        enodes = self.nodes_ee

        isize = enodes.shape[1] - 1
        jsize = enodes.shape[0]

        verts0 = onodes[:jsize-1, :isize]

        verts1 = enodes[1:jsize, :-1]
        verts2 = enodes[1:jsize, 1::]

        return _make_triangs(verts0, verts1, verts2)

    def triangs_from_verts2d(self, *verts):

        verts = [
            v.flatten() for v in verts
        ]

        return _pack_cols(*verts)


class _LatticeCloser:
    """Root of the lattice closer.
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self.cache = {}

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    @property
    def mesh_copy(self):
        return self.mesh.twin()

    @property
    def xsize(self):
        return self.mesh.meta['lattice-size'][0]

    @property
    def ysize(self):
        return self.mesh.meta['lattice-size'][1]

    @property
    def nodes0(self):
        return self.cache['side-nodes']['nodes-0']

    @property
    def nodes1(self):
        return self.cache['side-nodes']['nodes-1']

    @property
    def nodes3(self):
        return self.cache['side-nodes']['nodes-3']

    @property
    def nodes4(self):
        return self.cache['side-nodes']['nodes-4']

    @property
    def nodes2(self):
        return self.cache['new-points']['west']['nums']

    @property
    def nodes5(self):
        return self.cache['new-points']['east']['nums']

    EAST_534_SLICES = {
        (1, 0): (np.s_[:], np.s_[:], np.s_[:]),
        (1, 1): (np.s_[:], np.s_[:], np.s_[:-1]),
        (0, 1): (np.s_[1:], np.s_[1:], np.s_[:]),
        (0, 0): (np.s_[1:], np.s_[1:], np.s_[:-1])
    }

    EAST_354_SLICES = {
        (1, 0): (np.s_[:-1], np.s_[:-1], np.s_[1:]),
        (1, 1): (np.s_[:], np.s_[:], np.s_[1::]),
        (0, 1): (np.s_[:-1], np.s_[:-1], np.s_[:]),
        (0, 0): (np.s_[:], np.s_[:], np.s_[:])
    }


class LatticeCloser(_LatticeCloser):
    """Closes a primary lattice from the left.
    """

    def get_lattice_closed(self):

        self.make_new_points()
        self.make_new_triangs()

        mesh = self.mesh_copy

        mesh = mesh.add_points(
            self.cache['new-points']['west']['data']
        )

        mesh = mesh.add_points(
            self.cache['new-points']['east']['data']
        )

        mesh = mesh.add_triangs(
            self.cache['new-triangs']
        )

        mesh = mesh.renumed(
            np.argsort(mesh.points)
        )

        mesh = mesh.twin()
        mesh = mesh.add_meta(self.mesh.meta)

        return mesh

    def make_new_triangs(self):

        west = self.make_new_triangs_west()
        east = self.make_new_triangs_east()

        self.cache['new-triangs'] = np.vstack([west, east])

    def make_new_triangs_west(self):

        sz0 = self.nodes0.size - 1
        sz2 = self.nodes2.size

        triangs1 = _pack_cols(
            self.nodes0[:sz2], self.nodes1[:sz2], self.nodes2[:sz2]
        )

        triangs2 = _pack_cols(
            self.nodes2[:sz0], self.nodes1[:sz0], self.nodes0[1::]
        )

        return np.vstack([triangs1, triangs2])

    def make_new_triangs_east(self):

        key = (
            self.xsize % 2, self.ysize % 2
        )

        triangs1 = self.triangs_east_534(key)
        triangs2 = self.triangs_east_354(key)

        return np.vstack(
            [triangs1, triangs2]
        )

    def triangs_east_534(self, key):

        in3, in5, in4 = self.EAST_534_SLICES[key]

        return _pack_cols(
            self.nodes5[in5], self.nodes3[in3], self.nodes4[in4]
        )

    def triangs_east_354(self, key):

        in3, in5, in4 = self.EAST_354_SLICES[key]

        return _pack_cols(
            self.nodes3[in3], self.nodes5[in5], self.nodes4[in4]
        )

    def make_new_points(self):

        self.fetch_side_nodes()

        self.cache['new-points'] = {
            'west': self.make_new_points_west(),
            'east': self.make_new_points_east()
        }

    def make_new_points_west(self):

        data = self.mesh.points[self.nodes1] - 1.0
        bias = self.mesh.npoints

        return {
            'data': data,
            'nums': np.arange(data.size) + bias
        }

    def make_new_points_east(self):

        data = self.mesh.points[self.nodes3] + 1.0
        bias = self.mesh.npoints + self.nodes1.size

        return {
            'data': data,
            'nums': np.arange(data.size) + bias
        }

    def fetch_side_nodes(self):

        nodes0, = np.where(
            self.mesh.points.real == 0.0
        )

        nodes1, = np.where(
            self.mesh.points.real == 1.0
        )

        nodes3, = np.where(
            self.mesh.points.real == self.xsize - 2
        )

        nodes4, = np.where(
            self.mesh.points.real == self.xsize - 1
        )

        self.cache['side-nodes'] = {
            'nodes-0': nodes0,
            'nodes-1': nodes1,
            'nodes-3': nodes3,
            'nodes-4': nodes4
        }


def _pack_cols(*cols):
    return np.vstack(cols).T.copy('C')


def _make_triangs(*verts):

    verts = [
        v.flatten() for v in verts
    ]

    return _pack_cols(*verts)
