# -*- coding: utf-8 -*-
"""FEM unit partition.
"""
import numpy as np
from triellipt import geom


def getpartt(unit, loopnum, splitang, parttspec=None):
    """Creates a FEM unit partition.

    Parameters
    ----------
    unit : FEMUnit
        Unit to be partitioned.
    loopnum : int
        Number of the mesh loop to partition.
    splitang : float
        Threshold angle for splitting.
    parttspec : list = None
        Partition specification.

    Returns
    -------
    Partition
        Partition object.

    """

    partter = ParttMaker.from_unit(unit)

    partter = partter.with_params(
        loopnum, splitang, parttspec
    )

    return partter.get_partt()


class UnitAgent:
    """Operator on a FEM unit.
    """

    def __init__(self, unit):
        self.unit = unit
        self.cache = {}

    @classmethod
    def from_unit(cls, unit):
        return cls(unit)


class ParttMaker(UnitAgent):
    """Maker of a unit partition.
    """

    def with_params(self, loopnum, splitang, parttspec):

        self.cache['params'] = {
            'loopnum': loopnum,
            'splitang': splitang,
            'parttspec': parttspec
        }

        return self

    def get_partt(self):

        self.make_loop_colored()
        self.sync_loop_colored()

        paths = self.cache['loop-colored'].split()

        edges = {
            i + 1: v.numbers for i, v in enumerate(paths)
        }

        return FEMPartt(
            self.unit, edges
        )

    def sync_loop_colored(self):

        if not self.partition_spec:
            return

        coloroprs = {
            'l': self.lshift,
            'r': self.rshift,
            's': self.switch
        }

        for edgenum1, edgenum2, mode in self.partition_spec:
            coloroprs[mode](
                edgenum1 - 1, edgenum2 - 1
            )

    def make_loop_colored(self):

        cycle = self.make_loop_cycle()

        pathmap = cycle.split(
            self.cache['params']['splitang']
        )

        self.cache['loop-colored'] = pathmap

    def make_loop_cycle(self):

        cycle = geom.makecycle(
            self.take_target_loop()
        )

        self.cache['loop-cycle'] = cycle
        return cycle

    def take_target_loop(self):

        loop = self.unit.loops[
            self.cache['params']['loopnum']
        ]

        return loop.nodes_complex

    def lshift(self, color1, color2):
        self.pathmap = self.pathmap.lshift(color1, color2)

    def rshift(self, color1, color2):
        self.pathmap = self.pathmap.rshift(color1, color2)

    def switch(self, color1, color2):
        self.pathmap = self.pathmap.repaint(color1, color2)

    @property
    def pathmap(self):
        return self.cache['loop-colored']

    @pathmap.setter
    def pathmap(self, new_pathmap):
        self.cache['loop-colored'] = new_pathmap

    @property
    def partition_spec(self):
        return self.cache['params']['parttspec']


class FEMPartt:
    """FEM unit partition.

    Attributes
    ----------
    unit : FEMUnit
        Parent FEM unit.
    edge : dict
        Edge parts numbered from one.
    meta : dict
        Partition metadata.

    """

    def __init__(self, unit, edge):
        self.unit = unit
        self.edge = edge
        self.meta = self.push_meta()

    def push_meta(self):
        return {
            'dirich-edges': set()
        }

    @property
    def dirich_edges(self):
        return self.meta['dirich-edges']

    @property
    def neuman_edges(self):
        return set(self.edge) - self.dirich_edges

    def with_dirich_sides(self, *indices):
        self.meta['dirich-edges'] = set(indices)
        return self

    def edge2d(self, number):
        return self.unit.mesh.points2d[:, self.edge[number]]

    def asdict(self):
        return {
            0: self.core, **self.edge
        }

    @property
    def core(self):
        return self.getcore()

    def getcore(self):

        if not self.dirich_edges:
            return np.arange(self.unit.mesh_count)

        nodes_to_remove = [
            self.edge[ind] for ind in self.dirich_edges
        ]

        return np.delete(
            np.arange(self.unit.mesh_count), np.hstack(nodes_to_remove)
        )
