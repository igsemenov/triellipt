# -*- coding: utf-8 -*-
"""Stream of ij-coordinates.
"""
import numpy as np


def getstream(skeleton):
    """Returns the global matrix pattern.

    Returns
    -------
    IJStream
        Matrix pattern as an index stream.

    """
    return matrix_stream(skeleton)


def matrix_stream(skeleton):
    _ = MatrixStream.from_skeleton(skeleton)
    return _.get_stream()


def constr_stream(skeleton):
    _ = ConstrStream.from_skeleton(skeleton)
    return _.get_stream()


class SkelAgent:
    """Operator on a skeleton.
    """

    def __init__(self, skel):
        self.skel = skel

    @classmethod
    def from_skeleton(cls, skel):
        return cls(skel)

    @property
    def hasjoints(self):
        return self.skel.hasjoints

    @property
    def body(self):
        return self.skel.body

    @property
    def wests(self):
        return self.skel.wests

    @property
    def easts(self):
        return self.skel.easts

    @property
    def cores(self):
        return self.skel.cores

    @property
    def voids(self):
        return self.skel.voids


class MatrixStream(SkelAgent):
    """Matrix stream from a skeleton.
    """

    def get_stream(self):
        """Returns the full index stream.
        """

        body = self.get_stream_body()

        if not self.hasjoints:
            return body

        joints = self.get_stream_joints()

        return IJStream.from_data_train(body.data, joints.data)

    def get_stream_body(self):
        return self.get_stream_item('body')

    def get_stream_joints(self):

        wests = self.get_stream_wests()
        easts = self.get_stream_easts()
        cores = self.get_stream_cores()

        return IJStream.from_data_train(
            wests.data, easts.data, cores.data
        )

    def get_stream_wests(self):
        return self.get_stream_item('wests')

    def get_stream_easts(self):
        return self.get_stream_item('easts')

    def get_stream_cores(self):
        return self.get_stream_item('cores')

    def get_stream_item(self, key):

        offset = getattr(
            self, f'offset_{key}'
        )

        generator = getattr(
            self, f'gen_stream_{key}'
        )

        data = [
            stream.data for stream in generator(offset)
        ]

        return IJStream.from_data_train(*data)

    def gen_stream_body(self, offset):
        yield node_stream(self.body, self.body, 0, offset)
        yield node_stream(self.body, self.body, 1, offset)
        yield node_stream(self.body, self.body, 2, offset)

    def gen_stream_wests(self, offset):
        yield node_stream(self.wests, self.wests, 1, offset)
        yield node_stream(self.wests, self.wests, 1, offset)
        yield node_stream(self.wests, self.wests, 2, offset)

    def gen_stream_easts(self, offset):
        yield node_stream(self.easts, self.easts, 2, offset)
        yield node_stream(self.easts, self.easts, 1, offset)
        yield node_stream(self.easts, self.easts, 2, offset)

    def gen_stream_cores(self, offset):
        yield node_stream(self.wests, self.cores, 1, offset)
        yield node_stream(self.easts, self.cores, 2, offset)
        yield node_stream(self.cores, self.cores, 1, offset)
        yield node_stream(self.cores, self.cores, 2, offset)

    @property
    def offset_body(self):
        return 0

    @property
    def offset_wests(self):
        return self.body.size

    @property
    def offset_easts(self):
        return self.body.size + self.wests.size

    @property
    def offset_cores(self):
        return self.body.size + self.wests.size + self.easts.size


class ConstrStream(SkelAgent):
    """Constraint stream.
    """

    def get_stream(self):

        if not self.hasjoints:
            return None

        return node_stream(
            self.voids, self.voids, 2, self.offset_voids
        )

    @property
    def offset_voids(self):
        return self.body.size + sum(self.joints_sizes)

    @property
    def joints_sizes(self):
        return [
            self.wests.size, self.easts.size, self.cores.size
        ]


def node_stream(imesh, jmesh, node_index, trinums_offset):
    """Index stream from the specified local node.

    Parameters
    ----------
    imesh : TriMesh
        Mesh to take rows.
    jmesh : TriMesh
        Mesh to take cols.
    nodes_index : int
        Local index of the node in the rows-space.
    trinums_offset : int
        Offset of triangles numbers.

    """
    _ = NodeStream(imesh, jmesh)
    return _.get_stream(node_index, trinums_offset)


class IJStream:
    """Stream of ij-coordinates.
    """

    def __init__(self, data=None):
        self.data = data

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


class NodeStream:
    """Index stream from a node.

    Attributes
    ----------
    imesh : TriMesh
        Mesh to take rows.
    jmesh : TriMesh
        Mesh to take cols.

    """

    def __init__(self, imesh, jmesh):
        self.imesh = imesh
        self.jmesh = jmesh

    def get_stream(self, node_index, trinums_offset):

        cols = self.get_colnums_stream()
        rows = self.get_rownums_stream(node_index)
        tris = self.get_trinums_stream(trinums_offset)

        data = np.vstack(
            [rows, cols, tris]
        )

        return IJStream(data)

    def get_colnums_stream(self):
        return self.jmesh.triangs.flatten()

    def get_rownums_stream(self, node_index):
        return np.repeat(
            self.imesh.triangs[:, node_index], 3
        )

    def get_trinums_stream(self, trinums_offset):
        return np.repeat(self.trinums_range, 3) + trinums_offset

    @property
    def trinums_range(self):
        return np.arange(self.jmesh.size)
