# -*- coding: utf-8 -*-
"""FEM matrix factory.
"""
import numpy as np
from scipy import sparse as sp
from triellipt.utils import pairs


class FEMFactory:
    """Factory of FEM matrices.

    Attributes
    ----------
    mesh : TriMesh
        Parent mesh.
    body : sparse-matrix
        FEM-matrix pattern.
    meta : dict
        Factory metadata.

    """

    def __init__(self, unit=None, body=None, meta=None):
        self.unit = unit
        self.body = body
        self.meta = meta

    @classmethod
    def from_unit(cls, unit, add_constraints=True):
        return FactoryMaker(unit).get_factory(add_constraints)

    @property
    def shape(self):
        return self.body.shape

    @property
    def with_constraints(self):
        return self.meta['with-constraints']

    def __call__(self, data):
        return self.feed_data(data)

    def feed_data(self, data):
        """Transmits data to the FEM matrix.

        Parameters
        ----------
        data : flat-float-array
            Combination of local FEM operators (a).

        Returns
        -------
        csr-array
            Body of the resulting FEM matrix.
        dict
            Matrix metadata.

        Notes
        -----

        (a) Data stream compatible with ij-stream of the FEM unit.

        """

        body = self.make_body(data)
        meta = self.make_meta()

        return body, meta

    def make_body(self, data):

        data = self.push_data(data)

        new_body = self.body.copy()

        np.copyto(
            new_body.data[:self.frame_offset], data
        )

        return new_body

    def push_data(self, data):

        data_perm = self.meta['data-perm']
        bins_reduce = self.meta['bins-reduce']
        perm_reduced = self.meta['perm-reduced']

        data = np.add.reduceat(
            data[data_perm], bins_reduce
        )

        return np.copy(
            data[perm_reduced], order='C'
        )

    def make_meta(self):
        return {
            'has-constraints': self.with_constraints
        }

    @property
    def frame_offset(self):
        """Defines the length of actual data in the FEM matrix data.
        """
        return self.body.nnz - self.voids_offset

    @property
    def voids_offset(self):
        """Defines the length of constraints in the FEM matrix data.
        """
        if self.with_constraints:
            return 3 * self.unit.voids_count
        return 0

    @property
    def pattern(self):
        twin = self.body.copy()
        twin.data.fill(1)
        return twin


class UnitAgent:
    """Operator on a FEM unit.
    """

    def __init__(self, unit):
        self.unit = unit
        self.meta = self.fetch_meta()
        self.cache = self.push_cache()

    def fetch_meta(self):
        return {}

    def push_cache(self):
        return {}

    @classmethod
    def from_unit(cls, unit):
        return cls(unit)


class FactoryMaker(UnitAgent):
    """Maker of the FEM matrix factory.
    """

    def set_constraints_status(self, input_status):
        _ = self.get_constraints_status(input_status)
        self.cache['with-constraints'] = _

    def get_constraints_status(self, input_status):
        if not self.unit.hasvoids:
            return False
        return bool(
            input_status
        )

    @property
    def with_constraints(self):
        return self.cache['with-constraints']

    def get_factory(self, add_constraints):
        """Creates the FEM matrix factory.
        """

        self.set_constraints_status(add_constraints)

        ijmeta = self.make_ij_sorted_meta()
        matrix = self.from_ij_sorted_meta(ijmeta)

        return matrix

    def make_ij_sorted_meta(self):

        ij_sorter = IJSorter(
            *self.unit.ij_stream.ij_tuple
        )

        return ij_sorter.get_ij_sorted_meta()

    def from_ij_sorted_meta(self, ij_meta):

        body_data = self.make_body_data(ij_meta)
        meta_data = self.make_meta_data(ij_meta)

        return FEMFactory(
            self.unit, body_data, meta_data
        )

    def make_body_data(self, ij_meta):

        mold = self.make_matrix_mold(ij_meta)
        body = self.make_matrix_body(mold)

        return body

    def make_matrix_mold(self, ij_meta):

        ij_tuple = ij_meta['ij-tuple']

        index = np.arange(
            ij_tuple[0].size
        )

        mold_coo = sp.coo_array(
            (index, ij_tuple), shape=self.matrix_shape
        )

        mold_csr = mold_coo.tocsr()

        self.cache['perm-reduced'] = mold_csr.data.copy('C')
        return mold_csr

    def make_matrix_body(self, mold):

        body = mold.copy()

        body.data = np.ones_like(
            body.data, dtype=float
        )

        if not self.with_constraints:
            return body

        body = self.hang_constraint(body)
        return body

    def hang_constraint(self, body):
        if self.unit.hasvoids:
            return body + self.make_constraint()
        return body

    def make_constraint(self):
        _ = ConstrMaker.from_unit(self.unit)
        return _.get_matrix()

    def make_meta_data(self, ij_meta):

        data_perm = ij_meta['data-perm']
        bins_reduce = ij_meta['bins-reduce']

        perm_reduced = self.cache['perm-reduced']
        with_constraints = self.cache['with-constraints']

        return {
            'data-perm': data_perm,
            'bins-reduce': bins_reduce,
            'perm-reduced': perm_reduced,
            'with-constraints': with_constraints
        }

    @property
    def matrix_shape(self):
        return (
            self.unit.mesh_count, self.unit.mesh_count
        )


class Constrainer(UnitAgent):
    """Maker of constraints related matrices.
    """

    def fetch_meta(self):
        return {
            'voids-triangs': self.fetch_voids_triangs()
        }

    def fetch_voids_triangs(self):
        return self.unit.mesh.triangs[
            self.unit.mesh.meta['voids'], :
        ]

    def get_matrix(self):

        vij = self.make_v_ij()
        mat = self.from_v_ij(vij)

        return mat

    def make_ij_stream(self):

        i_stream = self.make_i_stream()
        j_stream = self.make_j_stream()

        return (i_stream, j_stream)

    def make_v_ij(self):

        v_stream = self.make_v_stream()
        ij_stream = self.make_ij_stream()

        return (
            v_stream, ij_stream
        )

    def from_v_ij(self, v_ij):

        mat_coo = sp.coo_array(
            v_ij, shape=self.matrix_shape
        )

        return mat_coo.tocsr()

    def make_i_stream(self):
        return ...

    def make_j_stream(self):
        return ...

    def make_v_stream(self):
        return ...

    @property
    def voids_triangs(self):
        return self.meta['voids-triangs']

    @property
    def matrix_shape(self):
        return (
            self.unit.mesh_count, self.unit.mesh_count
        )


class ConstrMaker(Constrainer):
    """Makes a constraint matrix.
    """

    def make_i_stream(self):
        return np.repeat(
            self.voids_triangs[:, 2], 3
        )

    def make_j_stream(self):
        return self.voids_triangs.flatten()

    def make_v_stream(self):
        return np.tile(
            self.row_constraint, self.unit.voids_count
        )

    @property
    def row_constraint(self):
        return np.array([1., 1., -2.])


class ProjerMaker(Constrainer):
    """Makes a projector to the constrained space.
    """

    def push_cache(self):
        return {
            'nodes-free': self.make_nodes_free(),
            'nodes-constr': self.make_nodes_constr()
        }

    def make_nodes_constr(self):
        return self.voids_triangs[:, 2]

    def make_nodes_free(self):
        return np.delete(
            np.arange(self.unit.mesh_count), self.voids_triangs[:, 2]
        )

    @property
    def nodes_free(self):
        return self.cache['nodes-free']

    @property
    def nodes_constr(self):
        return self.cache['nodes-constr']

    def make_i_stream(self):

        i_free = self.nodes_free

        i_constr = np.repeat(
            self.nodes_constr, 2
        )

        return np.hstack(
            [i_free, i_constr]
        )

    def make_j_stream(self):

        j_free = self.nodes_free
        j_constr = self.voids_triangs[:, [0, 1]].flatten()

        return np.hstack(
            [j_free, j_constr]
        )

    def make_v_stream(self):

        v_free = np.ones_like(self.nodes_free)

        v_constr = np.tile(
            [0.5, 0.5], self.unit.voids_count
        )

        return np.hstack(
            [v_free, v_constr]
        )


def get_constr_proj(unit):
    _ = ProjerMaker.from_unit(unit)
    return _.get_matrix()


class IJSorter:
    """Sorter of ij pairs.
    """

    def __init__(self, i_stream, j_stream):
        self.i_stream = i_stream
        self.j_stream = j_stream

    def get_ij_sorted_meta(self):

        data = self.make_ij_pairs_meta()
        meta = self.from_ij_pairs_meta(data)

        return meta

    def from_ij_pairs_meta(self, meta):

        data_perm = meta['pairs-sort']
        bins_reduce = meta['bins-reduce']
        packs_fronts = meta['packs-fronts']

        i_unique = self.i_stream[data_perm][packs_fronts]
        j_unique = self.j_stream[data_perm][packs_fronts]

        ij_tuple = (
            i_unique, j_unique
        )

        return {
            'ij-tuple': ij_tuple,
            'data-perm': data_perm,
            'bins-reduce': bins_reduce
        }

    def make_ij_pairs_meta(self):

        ij_pairs = self.make_ij_pairs()
        ij_metas = self.from_ij_pairs(ij_pairs)

        return ij_metas

    def make_ij_pairs(self):
        return pairs.szupaired(
            self.i_stream, self.j_stream
        )

    def from_ij_pairs(self, ij_pairs):
        """Returns ij-sorting metadata.
        """

        ij_sorted, sorter = _sort_data(ij_pairs)

        return {
            'pairs-sort': sorter, **_data_split(ij_sorted)
        }


def _sort_data(data):
    sorter = np.argsort(data).astype(int)
    return data[sorter], sorter


def _data_split(data):

    _, inds = np.unique(
        data, return_index=True
    )

    bins_split = inds[1::].astype(int)
    bins_reduce = np.hstack([0, bins_split])
    packs_fronts = np.hstack([0, bins_split])

    return {
        'bins-reduce': bins_reduce,
        'packs-fronts': packs_fronts
    }
