# -*- coding: utf-8 -*-
"""FEM unit partition.
"""
from abc import ABC, abstractmethod
import itertools as itr
import numpy as np
from triellipt import geom
from triellipt.fem import femvector, femmatrix


def getpartt(unit, spec):
    """Creates a FEM unit partition.

    Parameters
    ----------
    spec : dict
        Partition specification.

    Returns
    -------
    FEMPartt
        Partition object.

    """

    partt_name = spec['partition-title']
    loopspartt = spec['partition-loops']
    dirichlets = spec['dirichlet-sides']

    partt = makepartt(unit, loopspartt)
    partt = partt.with_name(partt_name)
    partt = partt.with_dirich_sides(*dirichlets)

    return partt


def makepartt(unit, loopspartt):
    _ = UnitPartter.from_unit(unit).with_loops_partt(loopspartt)
    return _.get_unit_partt()


class UnitAgent(ABC):
    """Operator on a FEM unit.
    """

    def __init__(self, unit):
        self.unit = unit
        self.cache = {}

    @classmethod
    def from_unit(cls, unit):
        return cls(unit)


class UnitPartter(UnitAgent):
    """Makes unit partition.
    """

    def with_loops_partt(self, loopspartt):
        self.cache['loops-partt'] = self.make_loops_partt(loopspartt)
        return self

    @property
    def loops_partt(self):
        return self.cache['loops-partt']

    def make_loops_partt(self, data):
        return LoopsPartt.from_attrs(self.unit, data).checked()

    def get_unit_partt(self):

        sides = self.get_sides_across_loops()

        sides = {
            i + 1: v for i, v in enumerate(sides)
        }

        return FEMPartt(
            self.unit, sides
        )

    def get_sides_across_loops(self):
        return itr.chain.from_iterable(
            self.gen_sides_across_loops()
        )

    def gen_sides_across_loops(self):

        for num, loop in enumerate(self.unit.loops):

            splitter = self.loops_partt.get_splitter_for_loop(num)
            coloring = self.loops_partt[num].get('coloring')

            yield splitter(loop, coloring)


LoopsParttError = type(
    'LoopsParttError', (Exception,), {}
)


class LoopsPartt:
    """Loops partition metadata.
    """

    def __init__(self, unit, data):
        self.unit = unit
        self.data = data

    @classmethod
    def from_attrs(cls, unit, data):
        return cls(unit, data)

    @property
    def loops_count(self):
        return len(self.unit.loops)

    def __getitem__(self, num):
        return self.data[num]

    def get_splitter_for_loop(self, loopnum):

        if 'angle' in self.data[loopnum]:
            return splitter_by_angle(
                self.data[loopnum]['angle']
            )

        if 'bins' in self.data[loopnum]:
            return splitter_by_bins(
                self.data[loopnum]['bins']
            )

        return None

    def checked(self):
        self.check_loops_numbers()
        self.check_partt_types()
        return self

    def check_loops_numbers(self):
        for loopnum in range(self.loops_count):
            if loopnum not in self.data:
                raise LoopsParttError(
                    f"no partition data for the loop {loopnum}"
                )

    def check_partt_types(self):
        for loopnum in range(self.loops_count):
            self.check_partt_type(loopnum)

    def check_partt_type(self, loopnum):
        if 'angle' in self.data[loopnum]:
            return
        if 'bins' in self.data[loopnum]:
            return
        raise LoopsParttError(
            f"no partition type for the loop {loopnum}"
        )


def splitter_by_bins(bins):
    def splitter_(loop, coloring):
        _ = SplitByBins(loop, coloring).with_bins(bins)
        return _.get_edges()
    return splitter_


def splitter_by_angle(angle):
    def splitter_(loop, coloring):
        _ = SplitByAngle(loop, coloring).with_angle(angle)
        return _.get_edges()
    return splitter_


class LoopSplitter:
    """Base splitter of loops.

    Parameters
    ----------
    loop : mesh-loop
        Loop to be partitioned.
    spec : dict
        Coloring specification.

    """

    def __init__(self, loop, spec):
        self.loop = loop
        self.spec = spec
        self.cache = {}

    def get_edges(self):

        self.make_loop_cycle()
        self.make_loop_colored()
        self.sync_loop_colored()

        paths = self.cache['loop-colored'].split()

        return [
            path.numbers for path in paths
        ]

    def sync_loop_colored(self):

        if not self.spec:
            return

        coloroprs = {
            'lshift': self.lshift,
            'rshift': self.rshift,
            'switch': self.switch
        }

        for sidenum1, sidenum2, mode in self.spec:
            coloroprs[mode](
                sidenum1 - 1, sidenum2 - 1
            )

    @abstractmethod
    def make_loop_colored(self):
        """Makes the pathmap from the loop-cycle.
        """

    def make_loop_cycle(self):

        cycle = geom.makecycle(
            self.loop.nodes_complex
        )

        self.cache['loop-cycle'] = cycle
        return cycle

    def lshift(self, color1, color2):
        self.pathmap = self.pathmap.lshift(color1, color2)

    def rshift(self, color1, color2):
        self.pathmap = self.pathmap.rshift(color1, color2)

    def switch(self, color1, color2):
        self.pathmap = self.pathmap.repaint(color1, color2)

    @property
    def cycle(self):
        return self.cache['loop-cycle']

    @property
    def pathmap(self):
        return self.cache['loop-colored']

    @pathmap.setter
    def pathmap(self, new_pathmap):
        self.cache['loop-colored'] = new_pathmap


class SplitByBins(LoopSplitter):
    """Splits the unit loop by bins.
    """

    def with_bins(self, bins):
        self.cache['split-bins'] = bins
        return self

    def make_loop_colored(self):
        pathmap = self.cycle.split(self.split_bins)
        self.cache['loop-colored'] = pathmap

    @property
    def split_bins(self):
        return self.cache['split-bins']


class SplitByAngle(LoopSplitter):
    """Splits the unit loop by an angle.
    """

    def with_angle(self, angle):
        self.cache['split-angle'] = angle
        return self

    def make_loop_colored(self):
        pathmap = self.cycle.dissect(self.split_angle)
        self.cache['loop-colored'] = pathmap

    @property
    def split_angle(self):
        return self.cache['split-angle']


def make_partt_base(unit):
    """Makes the base edge-core-partition.
    """

    edges = [
        loop.nodnums_unique for loop in unit.loops
    ]

    partt = FEMPartt(
        unit, {1: np.hstack(edges)}
    )

    partt = partt.with_name('base')
    partt = partt.with_dirich_sides(1)

    return partt


FEMParttError = type(
    'FEMParttError', (Exception,), {}
)


class DataPartt:
    """Partition data.

    Attributes
    ----------
    unit : FEMUnit
        Parent FEM unit.
    edge : dict
        Edge sides numbered from one.
    meta : dict
        Partition metadata.

    """

    def __init__(self, unit=None, edge=None, meta=None):
        self.unit = unit
        self.edge = edge
        self.meta = meta or {}

    def __getitem__(self, key):
        if key == 0:
            return self.core
        if key in self.edge:
            return self.edge[key]
        raise FEMParttError(
            f"no edge {key} in the partition '{self.name}'"
        )

    @property
    def name(self):
        if 'name' in self.meta:
            return self.meta['name']
        return str(self)

    @property
    def dirich_sides(self):
        return self.meta.get('dirich_sides')

    def update_meta(self, **kwargs):
        return self.__class__(
            self.unit, self.edge, self.meta | kwargs
        )

    @property
    def core(self):
        return self.getcore()

    def getcore(self):

        if not self.dirich_sides:
            return np.arange(
                self.unit.mesh_count
            )

        nodes_to_remove = [
            self.edge[ind] for ind in self.dirich_sides
        ]

        return np.delete(
            np.arange(self.unit.mesh_count), np.hstack(nodes_to_remove)
        )

    def nodes2d(self, key):
        return self.unit.mesh.points2d[:, self[key]]

    def nodes_complex(self, key):
        return self.unit.mesh.points[self[key]]


class FEMPartt(DataPartt):
    """FEM unit partition.

    Properties
    ----------

    Name        | Description
    ------------|----------------------------
    `name`      | Name of the partition.
    `edge`      | Map of the edge sections.
    `core`      | Core partition section.
    `meta`      | Partition metadata.

    """

    def with_name(self, name):
        return self.update_meta(name=name)

    def with_dirich_sides(self, *side_numbers):
        return self.update_meta(
            dirich_sides=set(side_numbers)
        )

    def new_vector(self, data=None):
        """Creates a new FEM vector.

        Parameters
        ----------
        data : scalar | flat-array
            Vector initialization data (optional).

        Returns
        -------
        VectorFEM

        """

        vector = femvector.getvector(self)

        if data is None:
            return vector
        return vector.with_body(data)

    def new_matrix(self, data, constr):
        """Creates a new FEM matrix.

        Parameters
        ----------
        data : flat-float-array
            Combination of local FEM operators.
        constr : bool
            Constraints are included, if True.

        Returns
        -------
        MatrixFEM
            Resulting FEM matrix.

        """
        if constr is False:
            return self.make_free_matrix(data)
        return self.make_full_matrix(data)

    def make_free_matrix(self, data):
        body, meta = self.unit.factory_free.feed_data(data)
        return femmatrix.getmatrix(self, body, meta)

    def make_full_matrix(self, data):
        body, meta = self.unit.factory_full.feed_data(data)
        return femmatrix.getmatrix(self, body, meta)

    def get_nodes(self, key):
        """Retrieves the points of the partition section.

        Parameters
        ----------
        key : int
            Number of the partition section.

        Returns
        -------
        two-row-float-array
            Points of the partition section stacked horizontally.

        """
        return self.unit.mesh.points2d[:, self[key]]
