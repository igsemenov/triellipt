# -*- coding: utf-8 -*-
"""Stream of matrix entries.
"""
from abc import abstractmethod, ABC
import numpy as np
from triellipt.fem import femoprs


def getstream(skeleton):
    """Returns a stream of FEM operators.
    """
    _ = FEMStream.from_skeleton(skeleton).with_operators()
    return _.get_fem_stream()


class SkelAgent:
    """Operator on a skeleton.
    """

    def __init__(self, skel):
        self.skel = skel
        self.oprs = None

    @classmethod
    def from_skeleton(cls, skel):
        return cls(skel)

    def with_operators(self):
        self.oprs = self.get_operators()
        return self

    def get_operators(self):

        body = {
            'body': femoprs.getoprs(self.skel.body)
        }

        if not self.skel.hasjoints:
            return body

        joints = {
            'wests': femoprs.getoprs(self.skel.wests),
            'easts': femoprs.getoprs(self.skel.easts),
            'cores': femoprs.getoprs(self.skel.cores)
        }

        return {
            **body, **joints
        }


class FEMStream(SkelAgent):
    """Maker of FEM streams.
    """

    def get_fem_stream(self):
        return self.get_matrix_stream()

    def get_matrix_stream(self):

        massmat = MassMatStream(self).get_stream()
        massdiag = MassDiagStream(self).get_stream()

        diff_1x = FluxStream(self, 'diff_1x').get_stream()
        diff_1y = FluxStream(self, 'diff_1y').get_stream()

        diff_2x = FluxStream(self, 'diff_2x').get_stream()
        diff_2y = FluxStream(self, 'diff_2y').get_stream()

        return {
            'massmat': massmat,
            'massdiag': massdiag,
            'diff_1x': diff_1x,
            'diff_1y': diff_1y,
            'diff_2x': diff_2x,
            'diff_2y': diff_2y
        }

    def get_constr_stream(self):

        if not self.skel.hasjoints:
            return None

        proxy = np.ones_like(
            self.skel.voids.triangs
        )

        return node_stream(
            proxy * [1., 1., -2.]
        )


class SubAgent(ABC):
    """Operator on a skeleton agent.
    """

    def __init__(self, agent):
        self.agent = agent

    @property
    def hasjoints(self):
        return self.agent.skel.hasjoints


class BaseStream(SubAgent):
    """Base operator streaming.
    """

    def get_stream(self):

        body = self.get_body()

        if not self.hasjoints:
            return body

        joints = self.get_joints()

        return np.hstack(
            [body, joints]
        )

    @abstractmethod
    def get_body(self):
        """Operator streaming from a body.
        """

    def get_joints(self):
        """Operator streaming from joints.
        """

        wests = self.get_wests()
        easts = self.get_easts()
        cores = self.get_cores()

        return np.hstack(
            [wests, easts, cores]
        )

    @abstractmethod
    def get_wests(self):
        """Operator streaming from wests.
        """

    @abstractmethod
    def get_easts(self):
        """Operator streaming from easts.
        """

    @abstractmethod
    def get_cores(self):
        """Operator streaming from cores.
        """

    @property
    def body(self):
        return self.agent.oprs['body']

    @property
    def wests(self):
        return self.agent.oprs['wests']

    @property
    def easts(self):
        return self.agent.oprs['easts']

    @property
    def cores(self):
        return self.agent.oprs['cores']


class MassMatStream(BaseStream):
    """Stream of the mass-matrix.
    """

    MASS_KEY = 'massmat'

    CORE_SCALE_WEST = np.array([0.5, 1.0, 0.0])
    CORE_SCALE_EAST = np.array([0.5, 0.0, 1.0])

    def get_body(self):
        return oprstream(
            self.body[self.MASS_KEY]
        )

    def get_wests(self):
        return oprstream(
            self.wests[self.MASS_KEY]
        )

    def get_easts(self):
        return oprstream(
            self.easts[self.MASS_KEY]
        )

    def get_cores(self):

        node0_west = node_stream(
            self.cores[self.MASS_KEY][0] * self.CORE_SCALE_WEST
        )

        node0_east = node_stream(
            self.cores[self.MASS_KEY][0] * self.CORE_SCALE_EAST
        )

        node1 = node_stream(self.cores[self.MASS_KEY][1])
        node2 = node_stream(self.cores[self.MASS_KEY][2])

        stream = [
            node0_west, node0_east, node1, node2
        ]

        return np.hstack(stream)


class MassDiagStream(MassMatStream):
    """Stream of the lumped mass-matrix.
    """

    MASS_KEY = 'massdiag'


class FluxStream(BaseStream):
    """Flux-operators streaming.
    """

    def __init__(self, agent, oprkey):
        super().__init__(agent)
        self.oprkey = oprkey

    def get_body(self):
        return oprstream(
            self.body[self.oprkey]
        )

    def get_wests(self):
        return oprstream(
            self.wests[self.oprkey]
        )

    def get_easts(self):
        return oprstream(
            self.easts[self.oprkey]
        )

    def get_cores(self):

        opr = self.cores[self.oprkey]

        flux0_west = -1. * opr[2]
        flux0_east = -1. * opr[1]

        node0_west = node_stream(flux0_west)
        node0_east = node_stream(flux0_east)

        node1 = node_stream(opr[1])
        node2 = node_stream(opr[2])

        stream = [
            node0_west, node0_east, node1, node2
        ]

        return np.hstack(stream)


def oprstream(operator):

    nodes = [
        node_stream(nodemat) for nodemat in operator
    ]

    return np.hstack(nodes)


def node_stream(node_matrix):
    return node_matrix.flatten()
