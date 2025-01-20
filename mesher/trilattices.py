# -*- coding: utf-8 -*-
"""Triangle lattice.
"""
import numpy as np
from triellipt.utils import pairs, tables
from triellipt.trimesh import trimesh_


def get_lattice(xsize, ysize):
    _ = get_latticer(xsize, ysize)
    return _.get_lattice()


def get_latticer(xsize, ysize):
    _ = TriLatticer()
    return _.with_counts(xsize, ysize)


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
            if self.icount >= self.jcount:
                return np.argsort(data)
            return np.argsort(
                np.conj(data * 1j)
            )
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
            'trilattice': {
                'xsize': self.icount,
                'ysize': self.jcount
            }
        }

        return mesh.add_meta(meta)

    def make_mesh(self):

        # Order matters: first triangles, then points.
        # Points renumbering depends on triangles.
        # Information transfer is via the cache.
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


def _pack_cols(*cols):
    return np.vstack(cols).T.copy('C')


def _make_triangs(*verts):

    verts = [
        v.flatten() for v in verts
    ]

    return _pack_cols(*verts)
