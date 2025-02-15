# -*- coding: utf-8 -*-
"""FEM matrix object.
"""


def getmatrix(partt, body, meta):
    """Creates a FEM matrix.
    """
    return MatrixFEM.from_data(partt, body, meta)


class MatrixData:
    """Root of the FEM matrix.
    """

    def __init__(self, partt=None, body=None, meta=None):
        self.partt = partt
        self.body = body
        self.meta = meta

    @classmethod
    def from_data(cls, unit, body, meta):
        return cls(unit, body, meta)

    @property
    def unit(self):
        return self.partt.unit

    @property
    def shape(self):
        return self.body.shape

    @property
    def has_constraints(self):
        return self.meta['has-constraints']

    def get_block_indexer(self, row_id, col_id):

        row_indexer = self.get_section_indexer(row_id)
        col_indexer = self.get_section_indexer(col_id)

        return row_indexer, col_indexer

    def get_section_indexer(self, key):
        return self.partt[key]

    def with_no_zeros(self):
        self.body.eliminate_zeros()
        return self

    @property
    def pattern(self):
        body = self.body_copy
        body.data.fill(1)
        return body

    @property
    def body_copy(self):
        return self.body.copy()


class MatrixFEM(MatrixData):
    """Global FEM matrix.
    """

    def __call__(self, row_key, col_key):
        return self.getblock(row_key, col_key)

    def __matmul__(self, vect):
        return vect.update_body(self.body @ vect.body)

    def getblock(self, row_id, col_id):
        """Extracts a block of a matrix.

        Parameters
        ----------
        row_id : int
            ID of the vertical section.
        col_id : int
            ID of the horizontal section.

        Returns
        -------
        csc-matrix
            Matrix bock in CSC format.

        """

        rows, cols = self.get_block_indexer(row_id, col_id)
        blockascsc = self.get_block_as_csc(rows, cols)

        return blockascsc

    def get_block_as_csc(self, rows, cols):

        body_ = self.body

        panel = body_[rows, :].tocsc()
        block = panel[:, cols]

        return block
