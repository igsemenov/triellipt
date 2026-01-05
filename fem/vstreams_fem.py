# -*- coding: utf-8 -*-
"""Stream of matrix entries (FEM).
"""
import numpy as np
from triellipt.fem import femoprs


def getstreams(skeleton):
    _ = StreamOprs.from_skel(skeleton)
    return _.get_streams()


class SkelAgent:
    """Operator on a skeleton.
    """

    def __init__(self, skel):
        self.skel = skel
        self.meta = self.fetch_meta()
        self.cache = {}

    @classmethod
    def from_skel(cls, skel):
        return cls(skel)

    def fetch_meta(self):
        return {}

    @property
    def hasvoids(self):
        return self.skel.hasvoids


class _StreamOprs(SkelAgent):
    """Collects streams of FEM opeators.
    """

    def fetch_meta(self):
        return {
            'fem-oprs': self.fetch_fem_oprs()
        }

    def fetch_fem_oprs(self):
        return femoprs.getoprs(self.skel.mesh)

    @property
    def massmat(self):
        return self.meta['fem-oprs']['massmat']

    @property
    def massdig(self):
        return self.meta['fem-oprs']['massdig']

    @property
    def diff_1x(self):
        return self.meta['fem-oprs']['diff_1x']

    @property
    def diff_1y(self):
        return self.meta['fem-oprs']['diff_1y']

    @property
    def grad_1y(self):
        return self.meta['fem-oprs']['grad_1y']

    @property
    def grad_1x(self):
        return self.meta['fem-oprs']['grad_1x']

    @property
    def diff_2x(self):
        return self.meta['fem-oprs']['diff_2x']

    @property
    def diff_2y(self):
        return self.meta['fem-oprs']['diff_2y']

    @property
    def diff_xy(self):
        return self.meta['fem-oprs']['diff_xy']

    @property
    def diff_yx(self):
        return self.meta['fem-oprs']['diff_yx']


class StreamOprs(_StreamOprs):
    """Collects streams of FEM opeators.
    """

    def get_streams(self):

        streams = {
            'massmat': self.get_stream_massmat(),
            'massdig': self.get_stream_massdig(),
            'diff_1x': self.get_stream_diff_1x(),
            'diff_1y': self.get_stream_diff_1y(),
            'diff_2x': self.get_stream_diff_2x(),
            'diff_2y': self.get_stream_diff_2y(),
            'diff_xy': self.get_stream_diff_xy(),
            'diff_yx': self.get_stream_diff_yx(),
            'grad_1y': self.get_stream_grad_1y(),
            'grad_1x': self.get_stream_grad_1x()
        }

        streams = {
            k: v.get_stream() for k, v, in streams.items()
        }

        return streams

    def get_stream_massmat(self):
        return self.stream_opr.with_opr(self.massmat)

    def get_stream_massdig(self):
        return self.stream_opr.with_opr(self.massdig)

    def get_stream_diff_1x(self):
        return self.stream_opr.with_opr(self.diff_1x)

    def get_stream_diff_1y(self):
        return self.stream_opr.with_opr(self.diff_1y)

    def get_stream_grad_1y(self):
        return self.stream_opr.with_opr(self.grad_1y)

    def get_stream_grad_1x(self):
        return self.stream_opr.with_opr(self.grad_1x)

    def get_stream_diff_2x(self):
        return self.stream_opr.with_opr(self.diff_2x)

    def get_stream_diff_2y(self):
        return self.stream_opr.with_opr(self.diff_2y)

    def get_stream_diff_xy(self):
        return self.stream_opr.with_opr(self.diff_xy)

    def get_stream_diff_yx(self):
        return self.stream_opr.with_opr(self.diff_yx)

    @property
    def stream_opr(self):
        return StreamOpr.from_skel(self.skel)


class StreamOpr(SkelAgent):
    """Base operator streamer.
    """

    def __init__(self, skel):
        super().__init__(skel)
        self.opr = None

    def with_opr(self, opr):
        self.opr = opr
        return self

    def get_stream(self):

        nodes_stream = self.get_stream_from_nodes()

        if not self.hasvoids:
            return nodes_stream

        voids_stream = self.get_stream_from_voids()

        return VStream.from_data_train(
            nodes_stream.data, voids_stream.data
        )

    def get_stream_from_nodes(self):
        return self.stream_from_map(self.nodesmap)

    def get_stream_from_voids(self):

        self.set_voids_submaps()

        west = self.stream_west()
        east = self.stream_east()

        stream = VStream.from_data_train(west.data, east.data)
        return stream

    def stream_west(self):

        west = self.stream_from_map(self.westmap)
        core = self.stream_from_map(self.coremap)
        east = self.stream_from_map(self.eastmap)

        return VStream.from_data_train(
            0.5 * west.data, 0.5 * core.data, 0.5 * east.data
        )

    def stream_east(self):

        west = self.stream_from_map(self.westmap)
        core = self.stream_from_map(self.coremap)
        east = self.stream_from_map(self.eastmap)

        return VStream.from_data_train(
            0.5 * west.data, 0.5 * core.data, 0.5 * east.data
        )

    def stream_from_map(self, srcmap):

        from_node_self = self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums + srcmap.locnums
        ]

        from_node_ccw1 = self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums + srcmap.locnums1
        ]

        from_node_ccw2 = self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums + srcmap.locnums2
        ]

        return VStream.from_data_train(
            from_node_self, from_node_ccw1, from_node_ccw2
        )

    def set_voids_submaps(self):
        self.cache |= self.skel.voids_submaps()

    @property
    def westmap(self):
        return self.cache['westmap']

    @property
    def coremap(self):
        return self.cache['coremap']

    @property
    def eastmap(self):
        return self.cache['eastmap']

    @property
    def nodesmap(self):
        return self.skel.nodesmap


class VStream:
    """Stream of matrix values.
    """

    def __init__(self, data):
        self.data = data

    @property
    def size(self):
        return self.data.size

    @classmethod
    def from_data(cls, data):
        return cls(data)

    @classmethod
    def from_data_train(cls, *data):
        return cls(np.hstack(data))

    def __mul__(self, value):
        return self.from_data(value * self.data)
