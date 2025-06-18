# -*- coding: utf-8 -*-
"""DtN unit.
"""
from triellipt import fem


def getdtn(mesh, anchors, mode="fem"):
    """Creates a DtN unit.

    Parameters
    ----------
    mesh : TriMesh
        Input triangle mesh.
    anchors : list[(float, float)]
        List of four anchor points.
    mode : str
        Solver mode — "fvm" or "fem" (default).

    Returns
    -------
    UnitDtN
        Resulting DtN unit.

    """
    return UnitDtN.from_mesh(mesh, anchors, mode)


class UnitDtN:
    """DtN unit.

    Attributes
    ----------
    unit : FEMUnit
        Background FEM unit.

    """

    def __init__(self, unit=None, mode=None):
        self.unit = unit
        self.mode = mode
        self.meta = {}

    @classmethod
    def from_mesh(cls, mesh, anchors, mode):
        return get_dtn_maker(mesh, anchors).make_dtn_unit(mode)

    @property
    def dtn(self):
        """DtN partition.
        """
        return self.unit.partts['dtn']

    @property
    def triu(self):
        return self.unit.mesh.triu

    @property
    def mesh(self):
        return self.unit.mesh

    @property
    def nodes1(self):
        return self.dtn.get_nodes(1)

    @property
    def nodes2(self):
        return self.dtn.get_nodes(2)

    @property
    def nodes3(self):
        return self.dtn.get_nodes(3)

    @property
    def nodes4(self):
        return self.dtn.get_nodes(4)

    @property
    def laplace(self):
        return self.dtn.new_matrix(
            self.unit.diff_2x + self.unit.diff_2y, add_constr=True
        )

    @property
    def massdiag(self):
        return self.dtn.new_matrix(
            self.unit.massdiag, add_constr=False
        )

    @property
    def massdiag_rad(self):
        return self.dtn.new_matrix(
            self.unit.radius * self.unit.massdiag, add_constr=False
        )

    @property
    def laplace_rad(self):

        lapl = self.unit.radius * (
            self.unit.diff_2x + self.unit.diff_2y
        )

        return self.dtn.new_matrix(lapl, add_constr=True)

    @property
    def laplace_axs(self):

        geom = self.unit.grad_1y / self.unit.radius

        return self.dtn.new_matrix(
            self.unit.diff_2x + self.unit.diff_2y + geom, add_constr=True
        )

    def new_vector(self):
        return self.dtn.new_vector()


def get_dtn_maker(mesh, anchors):
    return MakerDtN.from_mesh(mesh).with_anchors(anchors)


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self.unit = None
        self.meta = {}

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    def with_anchors(self, anchors: list):
        self.meta['anchors'] = anchors
        return self

    @property
    def anchors(self):
        return self.meta['anchors']


class MakerDtN(MeshAgent):
    """Basic DtN maker.
    """

    def make_dtn_unit(self, mode):

        self.with_fem_unit(mode)
        self.with_fem_unit_seeded()
        self.with_fem_unit_partted()

        return UnitDtN(self.unit, mode)

    def with_fem_unit(self, mode):

        self.unit = fem.getunit(
            self.mesh, anchors=[self.anchors[0]], mode=mode
        )

    def with_fem_unit_seeded(self):
        self.unit = self.unit.add_partition(self.partt_seed)

    def with_fem_unit_partted(self):
        self.unit = self.unit.add_partition(self.partt_dtn)

    def make_partt(self, name, anchor1, anchor2, anchor3):
        return {
            'name': name,
            'anchors': [
                anchor1, anchor2, anchor3
            ],
            'dirichlet-sides': [1, 2, 3, 4]
        }

    @property
    def partt_seed(self):
        return self.make_partt(
            'seed', self.anchors[1], self.anchors[2], self.anchors[3]
        )

    @property
    def partt_dtn(self):

        anchors = self.partt_dtn_anchors

        return self.make_partt(
            'dtn', anchors[0], anchors[1], anchors[2]
        )

    @property
    def partt_dtn_anchors(self):
        return self.get_dtn_anchors()

    def get_dtn_anchors(self):

        seed = self.unit.partts['seed']

        anchor1 = tuple(
            seed.get_nodes(2)[:, 1]
        )

        anchor2 = tuple(
            seed.get_nodes(3)[:, 0]
        )

        anchor3 = tuple(
            seed.get_nodes(4)[:, 1]
        )

        return [
            anchor1, anchor2, anchor3
        ]
