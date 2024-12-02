# -*- coding: utf-8 -*-
"""FEM matrix factory.
"""
import numpy as np
from scipy import sparse as sp
from triellipt.fem import femmatrix
from triellipt.utils import pairs


class FEMFactory:
    """Factory of FEM matrices.
    """

    def __init__(self, unit=None, body=None, meta=None):
        self.unit = unit
        self.body = body
        self.meta = meta

    @classmethod
    def from_unit(cls, unit, with_constraints=True):
        return FactoryMaker(unit).get_factory(with_constraints)

    @property
    def shape(self):
        return self.body.shape

    @property
    def frame_offset(self):
        return self.body.nnz - self.voids_offset

    @property
    def voids_offset(self):
        if self.has_constraints:
            return 3 * self.unit.voids_count
        return 0

    @property
    def has_constraints(self):
        return self.meta['has-constraints']

    def feed_data(self, data):
        """Transmits data to the matrix.

        Parameters
        ----------
        data : flat-float-array
            Data stream compatible with ij-stream of the FEM unit.

        Returns
        -------
        MatrixFEM
            Resulting FEM matrix.

        """

        new_data = self.push_data(data)
        fem_matr = self.make_matrix(new_data)

        return fem_matr

    def push_data(self, data):

        data = np.add.reduceat(
            data[self.data_perm], self.bins_reduce
        )

        return np.copy(
            data[self.perm_reduced], order='C'
        )

    def make_matrix(self, new_data):

        body = self.make_body(new_data)
        meta = self.make_meta()

        return femmatrix.getmatrix(self.unit, body, meta)

    def make_body(self, new_data):

        new_body = self.body.copy()

        np.copyto(
            new_body.data[:self.frame_offset], new_data
        )

        return new_body

    def make_meta(self):
        return {
            'has-constraints': self.has_constraints
        }

    def pattern(self):

        twin = self.body.copy()
        twin.data.fill(1)

        return twin

    @property
    def data_perm(self):
        """Input data permutation.
        """
        return self.meta['data-perm']

    @property
    def bins_reduce(self):
        """Bins to reduce input data.
        """
        return self.meta['bins-reduce']

    @property
    def perm_reduced(self):
        """Permutation of the reduced data.
        """
        return self.meta['perm-reduced']


class UnitAgent:
    """Operator on a FEM unit.
    """

    def __init__(self, unit):
        self.unit = unit

    @classmethod
    def from_unit(cls, unit):
        return cls(unit)

    @property
    def hasjoints(self):
        return self.unit.hasjoints


class FactoryMaker(UnitAgent):
    """Maker of the FEM matrix factory.
    """

    def __init__(self, unit):
        super().__init__(unit)
        self.with_constraints = None

    def get_factory(self, with_constraints):
        self.with_constraints = with_constraints

        ijmeta = self.make_ij_sorted_meta()
        matrix = self.from_ij_sorted_meta(ijmeta)

        return matrix

    def make_ij_sorted_meta(self):
        return self.ij_sorter.get_ij_sorted_meta()

    def from_ij_sorted_meta(self, ij_meta):

        body_data = self.make_body_data(ij_meta)
        meta_data = self.make_meta_data(ij_meta, body_data)

        return FEMFactory(
            self.unit, body_data['body'], meta_data
        )

    def make_body_data(self, ij_meta):

        mold = self.make_matrix_mold(ij_meta)
        body = self.make_matrix_body(mold)

        perm_reduced = mold.data.copy('C')

        return {
            'body': body, 'perm-reduced': perm_reduced
        }

    def make_meta_data(self, ij_meta, body_data):

        with_constraints = self.with_constraints

        return {
            'data-perm': ij_meta['data-perm'],
            'bins-reduce': ij_meta['bins-reduce'],
            'perm-reduced': body_data['perm-reduced'],
            'has-constraints': with_constraints
        }

    def make_matrix_mold(self, ij_meta):

        ij_tuple = ij_meta['ij-tuple']

        index = np.arange(
            ij_tuple[0].size
        )

        mold_coo = sp.coo_array(
            (index, ij_tuple), shape=self.matrix_shape
        )

        mold_csr = mold_coo.tocsr()
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
        if self.hasjoints:
            return body + self.make_constraint()
        return body

    def make_constraint(self):
        _ = ConstrMaker.from_unit(self.unit)
        return _.get_constr()

    @property
    def ij_sorter(self):
        return IJSorter(
            *self.unit.ij_stream.ij_tuple
        )

    @property
    def matrix_shape(self):
        return (
            self.unit.mesh_count, self.unit.mesh_count
        )


class ConstrMaker(UnitAgent):
    """Makes a constraint matrix.
    """

    def get_constr(self):

        vij = self.make_v_ij()
        mat = self.from_v_ij(vij)

        return mat

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

    def make_ij_stream(self):

        voids_triangs = self.voids_triangs()

        i_stream = self.make_i_stream(voids_triangs)
        j_stream = self.make_j_stream(voids_triangs)

        return (i_stream, j_stream)

    def make_i_stream(self, voids_triangs):
        return np.repeat(
            voids_triangs[:, 2], 3
        )

    def make_j_stream(self, voids_triangs):
        return voids_triangs.flatten()

    def make_v_stream(self):
        return np.tile(
            self.row_constraint, self.unit.voids_count
        )

    def voids_triangs(self):
        return self.unit.mesh_voids.triangs

    @property
    def row_constraint(self):
        return np.array([1., 1., -2.])

    @property
    def matrix_shape(self):
        return (
            self.unit.mesh_count, self.unit.mesh_count
        )


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

    def make_ij_pairs_meta(self):

        ij_pairs = self.make_ij_pairs()
        ij_metas = self.from_ij_pairs(ij_pairs)

        return ij_metas

    def from_ij_pairs_meta(self, meta):

        data_perm = meta['data-perm']
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

    def make_ij_pairs(self):
        return pairs.szupaired(
            self.i_stream, self.j_stream
        )

    def from_ij_pairs(self, ij_pairs):
        """Returns ij-sorting metadata.
        """

        data_perm = _argsort(ij_pairs)

        _, bins = np.unique(
            ij_pairs[data_perm], return_index=True
        )

        bins_split = bins[1::].astype(int)
        bins_reduce = np.r_[0, bins_split]
        packs_fronts = np.r_[0, bins_split]

        return {
            'data-perm': data_perm,
            'bins-reduce': bins_reduce,
            'packs-fronts': packs_fronts
        }


def _argsort(data):
    return np.argsort(data).astype(int)
