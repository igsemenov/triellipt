# -*- coding: utf-8 -*-
"""Stream of matrix entries (FVM).
"""
from abc import abstractmethod, ABC
import numpy as np
from triellipt.fem import femoprs


def getstreams(skeleton):
    _ = StreamOprs.from_skel(skeleton)
    return _.get_streams()


class SkelAgent(ABC):
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


class StreamOprs(SkelAgent):
    """Collects streams of FEM opeators.
    """

    def fetch_meta(self):
        return {
            'fem-oprs': self.fetch_fem_oprs()
        }

    def fetch_fem_oprs(self):
        return femoprs.getoprs(self.skel.mesh)

    def get_streams(self):

        streams = {
            'massmat': self.get_stream_massmat(),
            'massdiag': self.get_stream_massdiag(),
            'diff_1x': self.get_stream_diff_1x(),
            'diff_1y': self.get_stream_diff_1y(),
            'diff_2x': self.get_stream_diff_2x(),
            'diff_2y': self.get_stream_diff_2y(),
            'grad_1y': self.get_stream_grad_1y()
        }

        streams = {
            k: v.get_stream() for k, v, in streams.items()
        }

        return streams

    def get_stream_massmat(self):
        return self.stream_mass.with_opr(self.massmat)

    def get_stream_massdiag(self):
        return self.stream_mass.with_opr(self.massdiag)

    def get_stream_diff_1x(self):
        return self.stream_flux.with_opr(self.diff_1x)

    def get_stream_diff_1y(self):
        return self.stream_flux.with_opr(self.diff_1y)

    def get_stream_grad_1y(self):
        return self.stream_mass.with_opr(self.grad_1y)

    def get_stream_diff_2x(self):
        return self.stream_flux.with_opr(self.diff_2x)

    def get_stream_diff_2y(self):
        return self.stream_flux.with_opr(self.diff_2y)

    @property
    def stream_mass(self):
        return StreamMass.from_skel(self.skel)

    @property
    def stream_flux(self):
        return StreamFlux.from_skel(self.skel)

    @property
    def massmat(self):
        return self.meta['fem-oprs']['massmat']

    @property
    def massdiag(self):
        return self.meta['fem-oprs']['massdiag']

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
    def diff_2x(self):
        return self.meta['fem-oprs']['diff_2x']

    @property
    def diff_2y(self):
        return self.meta['fem-oprs']['diff_2y']


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

    def get_stream_from_voids(self):

        self.set_voids_submaps()

        west_west = self.stream_west_west()
        core_west = self.stream_core_west()
        core_east = self.stream_core_east()
        east_east = self.stream_east_east()

        stream = VStream.from_data_train(
            west_west.data,
            core_west.data,
            core_east.data,
            east_east.data
        )

        return stream

    def set_voids_submaps(self):
        self.cache |= self.skel.voids_submaps()

    def stream_west_west(self):
        return self.stream_from_map(self.westmap)

    def stream_east_east(self):
        return self.stream_from_map(self.eastmap)

    @abstractmethod
    def stream_core_west(self):
        """Stream from the core-west splitting.
        """

    @abstractmethod
    def stream_core_east(self):
        """Stream from the core-east splitting.
        """

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


class StreamMass(StreamOpr):
    """Mass operator streamer.
    """

    def stream_core_west(self):

        srcmap = self.coremap

        from_node_self = self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums + srcmap.locnums
        ]

        from_node_ccw1 = self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums + srcmap.locnums1
        ]

        from_node_ccw2 = np.zeros_like(from_node_self)

        return VStream.from_data_train(
            0.5 * from_node_self, from_node_ccw1, from_node_ccw2
        )

    def stream_core_east(self):

        srcmap = self.coremap

        from_node_self = self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums + srcmap.locnums
        ]

        from_node_ccw2 = self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums + srcmap.locnums2
        ]

        from_node_ccw1 = np.zeros_like(from_node_self)

        return VStream.from_data_train(
            0.5 * from_node_self, from_node_ccw1, from_node_ccw2
        )


class StreamFlux(StreamOpr):
    """Flux operator streamer.
    """

    def stream_core_west(self):

        srcmap = self.coremap

        from_node_self = -1. * self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums2 + srcmap.locnums
        ]

        from_node_ccw1 = -1. * self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums2 + srcmap.locnums1
        ]

        from_node_ccw2 = -1. * self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums2 + srcmap.locnums2
        ]

        return VStream.from_data_train(
            from_node_self, from_node_ccw1, from_node_ccw2
        )

    def stream_core_east(self):

        srcmap = self.coremap

        from_node_self = -1. * self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums1 + srcmap.locnums
        ]

        from_node_ccw1 = -1. * self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums1 + srcmap.locnums1
        ]

        from_node_ccw2 = -1. * self.opr.flat[
            9 * srcmap.trinums + 3 * srcmap.locnums1 + srcmap.locnums2
        ]

        return VStream.from_data_train(
            from_node_self, from_node_ccw1, from_node_ccw2
        )


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
