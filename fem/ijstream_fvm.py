# -*- coding: utf-8 -*-
"""Stream of ij-coordinates (FVM).
"""
import numpy as np


def getstream(skeleton):
    """Stream maker called by the FEM unit.
    """
    return matrix_stream(skeleton)


def matrix_stream(skeleton):
    """Creates the FEM matrix pattern.

    Parameters
    ----------
    skeleton : MeshSkeleton
        Mesh in skeleton format.

    Returns
    -------
    IJStream
        Matrix pattern as an ij-stream.

    """
    _ = MatrixStream.from_skeleton(skeleton)
    return _.get_stream()


class SkelAgent:
    """Operator on a mesh skeleton.
    """

    def __init__(self, skel):
        self.skel = skel
        self.cache = {}

    @classmethod
    def from_skeleton(cls, skel):
        return cls(skel)

    @property
    def hasvoids(self):
        return self.skel.hasvoids


class MatrixStream(SkelAgent):
    """Matrix stream from a skeleton.
    """

    def get_stream(self):

        nodes_stream = self.get_stream_from_nodes()

        if not self.hasvoids:
            return nodes_stream

        voids_stream = self.get_stream_from_voids()

        stream = IJStream.from_data_train(
            nodes_stream.data, voids_stream.data
        )

        stream.meta = {
            **nodes_stream.meta,
            **voids_stream.meta
        }

        return stream

    def get_stream_from_nodes(self):

        stream = self.stream_from_map(
            self.nodesmap.nodnums, self.nodesmap
        )

        stream.meta = {
            'hasvoids': False,
            'nodsmap-size': self.nodesmap.size
        }

        return stream

    def get_stream_from_voids(self):

        self.set_voids_submaps()

        west_west = self.stream_west_west()
        core_west = self.stream_core_west()
        core_east = self.stream_core_east()
        east_east = self.stream_east_east()

        stream = IJStream.from_data_train(
            west_west.data,
            core_west.data,
            core_east.data,
            east_east.data
        )

        stream.meta = {
            'hasvoids': True,
            'westmap-size': self.westmap.size,
            'coremap-size': self.coremap.size,
            'eastmap-size': self.eastmap.size
        }

        return stream

    def set_voids_submaps(self):
        self.cache |= self.skel.voids_submaps()

    def stream_west_west(self):
        return self.stream_from_map(
            self.westmap.nodnums1, self.westmap
        )

    def stream_core_west(self):
        return self.stream_from_map(
            self.westmap.nodnums1, self.coremap
        )

    def stream_core_east(self):
        return self.stream_from_map(
            self.eastmap.nodnums2, self.coremap
        )

    def stream_east_east(self):
        return self.stream_from_map(
            self.eastmap.nodnums2, self.eastmap
        )

    @classmethod
    def stream_from_map(cls, dstnodes, srcmap):

        from_node_self = _stream_node(
            dstnodes, srcmap.nodnums, srcmap.trinums
        )

        from_node_ccw1 = _stream_node(
            dstnodes, srcmap.nodnums1, srcmap.trinums
        )

        from_node_ccw2 = _stream_node(
            dstnodes, srcmap.nodnums2, srcmap.trinums
        )

        return IJStream.from_data_train(
            from_node_self, from_node_ccw1, from_node_ccw2
        )

    @property
    def nodesmap(self):
        return self.skel.nodesmap

    @property
    def westmap(self):
        return self.cache['westmap']

    @property
    def coremap(self):
        return self.cache['coremap']

    @property
    def eastmap(self):
        return self.cache['eastmap']


def _stream_node(dstnodes, srcnodes, trinums):

    rownums = dstnodes
    colnums = srcnodes

    return _pack_rows(
        rownums, colnums, trinums
    )


class IJStream:
    """Stream of ij-coordinates.
    """

    def __init__(self, data=None):
        self.data = data
        self.meta = {}

    @classmethod
    def from_data(cls, data):
        return cls(data)

    @classmethod
    def from_data_train(cls, *data):
        return cls.from_data(np.hstack(data))

    @property
    def size(self):
        return self.data.shape[1]

    @property
    def hasvoids(self):
        return self.meta['hasvoids']

    @property
    def rownums(self):
        return self.data[0, :]

    @property
    def colnums(self):
        return self.data[1, :]

    @property
    def trinums(self):
        return self.data[2, :]

    @property
    def ij_tuple(self):
        return (
            self.rownums, self.colnums
        )

    def update_data(self, new_data):
        return self.from_data(
            np.copy(new_data, order='C')
        )

    def atrow(self, rownum):
        """Returns a stream from the specified row.
        """

        new_data = self.data[
            :, self.rownums == rownum
        ]

        return self.update_data(new_data)

    def atcol(self, colnum):
        """Returns a stream from the specified col.
        """

        new_data = self.data[
            :, self.colnums == colnum
        ]

        return self.update_data(new_data)

    def mask_at_row(self, rownum):
        return self.rownums == rownum

    def mask_at_col(self, colnum):
        return self.colnums == colnum


def _pack_rows(*rows):
    return np.vstack(rows).copy('C')


def _pack_cols(*cols):
    return np.vstack(cols).T.copy('C')
