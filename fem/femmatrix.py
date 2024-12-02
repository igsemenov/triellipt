# -*- coding: utf-8 -*-
"""FEM matrix object.
"""
import numpy as np


def getmatrix(unit, body, meta):
    """Creates a FEM matrix.
    """
    return MatrixFEM.from_attrs(unit, body, meta)


class MatrixData:
    """Root of the FEM matrix.
    """

    def __init__(self, unit=None, body=None, meta=None):
        self.unit = unit
        self.body = body
        self.meta = meta

    @classmethod
    def from_attrs(cls, unit, body, meta):
        return cls(unit, body, meta)

    @property
    def shape(self):
        return self.body.shape

    @property
    def is_named(self):
        return 'name' in self.meta

    @property
    def haspartition(self):
        return 'partition' in self.meta

    @property
    def hasconstraints(self):
        return self.meta['has-constraints']

    @property
    def partition_name(self):
        return self.meta['partition']['name']

    @property
    def matrix_name(self):
        if self.is_named:
            return self.meta['name']
        return str(id(self))

    def update_meta(self, new_meta):
        return self.__class__(
            self.unit, self.body_copy, new_meta
        )

    def get_block_indexer(self, row_key, col_key):

        row_indexer = self.get_section_indexer(row_key)
        col_indexer = self.get_section_indexer(col_key)

        return row_indexer, col_indexer

    def get_section_indexer(self, name):

        if not self.haspartition:
            raise MatrixFEMError(
                f"matrix '{self.matrix_name}' has not partition"
            )

        indexer = self.meta['partition']['sections'].get(name)

        if indexer is None:
            raise MatrixFEMError(
                f"no section '{name}' in matrix '{self.matrix_name}'"
            )

        return indexer

    @property
    def body_copy(self):
        return self.body.copy()


class MatrixFEM(MatrixData):
    """Global FEM matrix.
    """

    def __matmul__(self, vect):
        """Multiplication by a FEM vector only.
        """
        return vect.update_body(
            self.body @ vect.body
        )

    def with_name(self, name):
        """Prescribes the matrix name.

        Parameters
        ----------
        name : str
            Name of the matrix.

        Returns
        -------
        MatrixFEM
            Copy of the matrix with the new name.

        """
        return self.add_meta({'name': name})

    def dirichsplit(self):
        """Creates a Dirichlet (core-edge) partition.
        """
        return self.partitioned(
            meta_edge_core(self.unit)
        )

    def partitioned(self, meta):
        """Creates a partitioned matrix.

        Parameters
        ----------
        meta : dict
            Partition meta data (name and sections).

        Returns
        -------
        MatrixFEM
            Copy of the matrix with the partition defined.

        """

        new_meta = {
            'partition': meta
        }

        return self.add_meta(new_meta)

    def add_meta(self, meta):
        return self.update_meta(self.meta | meta)

    def getblock(self, rowname, colname):
        """Extracts a block from a partitioned matrix.

        Parameters
        ----------
        rowname : str
            Name of the vertical section.
        colname : str
            Name of the horizontal section.

        Returns
        -------
        csc-matrix
            Block of a partitioned matrix.

        """

        rows, cols = self.get_block_indexer(rowname, colname)
        blockascsc = self.get_block_as_csc(rows, cols)

        return blockascsc

    def get_block_as_csc(self, rows, cols):

        body_ = self.body

        panel = body_[rows, :].tocsc()
        block = panel[:, cols]

        return block


MatrixFEMError = type(
    'MatrixFEMError', (Exception,), {}
)


def meta_edge_core(unit):
    """Dirichlet partition meta data for a FEM unit.
    """

    mesh_count = unit.mesh_count
    edge_count = unit.edge_count

    mesh_range = np.arange(mesh_count)

    edge_range = mesh_range[:edge_count]
    core_range = mesh_range[edge_count:]

    return {
        'name': 'edge-core',
        'sections': {
            'edge': edge_range,
            'core': core_range
        }
    }
