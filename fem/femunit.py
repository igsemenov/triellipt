# -*- coding: utf-8 -*-
"""Public FEM unit.
"""
import numpy as np
from triellipt.fem import (
    skeleton,
    ijstream,
    vstreams,
    femfactory,
    femoprs,
    fempartt,
    massinv_,
    trinterp
)


def getunit(mesh, anchors=None, mode=None):
    """Creates a FEM computing unit.

    Parameters
    ----------
    mesh : TriMesh
        Input triangle mesh.
    anchors : Iterable = None
        Provides `(float, float)` points to synchronize the mesh boundary.
    mode : str = None
        Solver mode — "fvm" or "fem" (default).

    Returns
    -------
    FEMUnit
        FEM computing unit.

    """
    return FEMUnit.from_mesh(
        mesh, anchors, mode or "fem"
    )


class FEMData:
    """Base FEM data.
    """

    def __init__(self, mesh=None, meta=None):
        self.mesh = mesh
        self.meta = meta
        self.cache = {}

    @classmethod
    def from_mesh(cls, mesh, anchors=None, mode="fem"):

        unit_data = FEMUnitMaker().get_unit_data(mesh, anchors, mode)

        unit = cls(
            unit_data["mesh"], unit_data["meta"]
        )

        return unit.with_base_partition()

    @property
    def hasvoids(self):
        return self.voids_count != 0

    @property
    def mesh_count(self):
        """Number of mesh nodes.
        """
        return self.mesh.npoints

    @property
    def voids_count(self):
        """Number of voids nodes.
        """
        return self.mesh.meta['voids'].size

    @property
    def frame_count(self):
        """Number of frame nodes.
        """
        return self.mesh_count - self.voids_count

    @property
    def edge_count(self):
        """Number of edge nodes.
        """

        sizes = [
            loop.size for loop in self.loops
        ]

        return sum(sizes)

    @property
    def core_count(self):
        """Number of core nodes.
        """
        return self.mesh_count - self.edge_count

    @property
    def loops(self):
        """Loops of the unit mesh.
        """
        return self.mesh.meta['loops']

    def with_base_partition(self):
        """Adds the default partition to the unit.
        """

        self.cache['partitions'] = {
            'base': fempartt.make_base_partition(self)
        }

        return self

    @property
    def femoprs(self):
        return self.meta['v-streams']

    @property
    def ij_stream(self):
        return self.meta['ij-stream']

    @property
    def ij_tuple(self):
        return (
            self.ij_stream.rownums, self.ij_stream.colnums
        )

    @property
    def ij_t(self):
        return self.ij_stream.trinums


class FEMRoot(FEMData):
    """Root of the FEM unit.
    """

    @property
    def massmat(self):
        return self.femoprs['massmat'].data

    @property
    def massdiag(self):
        return self.femoprs['massdiag'].data

    @property
    def diff_1x(self):
        return - self.femoprs['diff_1x'].data

    @property
    def diff_1y(self):
        return - self.femoprs['diff_1y'].data

    @property
    def grad_1y(self):
        return + self.femoprs['grad_1y'].data

    @property
    def diff_2x(self):
        return - self.femoprs['diff_2x'].data

    @property
    def diff_2y(self):
        return - self.femoprs['diff_2y'].data

    @property
    def radius(self):
        return self.mesh.centrs_complex.imag[self.ij_t]

    @property
    def grad(self):
        """Gradient operator.
        """

        if 'grad' in self.cache:
            return self.cache['grad']

        self.cache['grad'] = femoprs.getgrad(self.mesh)
        return self.grad

    @property
    def trigeo(self):
        """Dictionary with the geometric properties of triangles.
        """

        if 'trigeo' in self.cache:
            return self.cache['trigeo']

        self.cache['trigeo'] = femoprs.mesh_geom(self.mesh)
        return self.trigeo

    @property
    def perm(self):
        """Nodes permutation in the mesh-to-unit transition.
        """
        return self.meta['mesh2unit']

    @property
    def base(self):
        return self.partts['base']

    @property
    def partts(self):
        return self.cache['partitions']

    def fem_factory(self, add_constr=True):
        """Creates a factory of FEM matrices.

        Parameters
        ----------
        add_constr : bool = True
            If True, forces constraints to be included in the matrix.

        Returns
        -------
        FEMFactory
            Callable factory of FEM matrices.

        """
        return femfactory.FEMFactory.from_unit(self, add_constr)

    def massinv(self):
        """Creates the inverse mass operator.

        Returns
        -------
        MassDiagInv
            Callable inverse operator.

        Notes
        -----

        Used only for a lumped mass operator with constraints.

        """
        return massinv_.getmassinv(self)

    @property
    def factory_free(self):
        """Factory with no constraints.
        """

        if 'factory-free' in self.cache:
            return self.cache['factory-free']

        _ = self.fem_factory(add_constr=False)
        self.cache['factory-free'] = _
        return self.factory_free

    @property
    def factory_full(self):
        """Factory with constraints.
        """

        if 'factory-full' in self.cache:
            return self.cache['factory-full']

        _ = self.fem_factory(add_constr=True)
        self.cache['factory-full'] = _
        return self.factory_full


class FEMUnit(FEMRoot):
    """FEM computing unit.

    Properties
    ----------

    FEM operators as data-streams:

    Name        | Description
    ------------|----------------------
    `massmat`   | Mass-matrix
    `massdiag`  | Mass-matrix lumped
    `diff_1y`   | 1st y-derivative
    `diff_1x`   | 1st x-derivative
    `diff_2y`   | 2nd y-derivative
    `diff_2x`   | 2nd x-derivative

    General properties:

    Name      | Description
    ----------|-------------------------------
    `grad`    | Gradient operator.
    `trigeo`  | Geometric properties.
    `perm`    | Mesh-to-unit permutation.
    `base`    | Base edge-core partition.
    `loops`   | List of the mesh loops.
    `partts`  | Map of the unit partitions.

    """

    def add_partition(self, partt_spec):
        """Adds new partition to the unit.

        Parameters
        ----------
        partt_spec : dict
            Partition specification.

        Returns
        -------
        self
            Unit with the partition added.

        """

        self.set_partition(
            fempartt.getpartt(self, partt_spec)
        )

        return self

    def get_partition(self, partt_name):
        """Fetches the unit partition.

        Parameters
        ----------
        partt_name : str
            Partition name.

        Returns
        -------
        FEMPartt
            Desired unit partition.

        """
        if partt_name not in self.partts:
            raise ValueError(f"unit has no '{partt_name}' partition")
        return self.partts[partt_name]

    def set_partition(self, partt) -> None:
        """Assigns the partition to the unit.

        Parameters
        ----------
        partt : FEMPartt
            Input unit partition.

        """

        if partt.unit is not self:
            raise ValueError(
                f"partition '{partt.name}' is from another unit"
            )

        self.cache['partitions'][partt.name] = partt

    def del_partition(self, name) -> None:
        """Deletes the specified partition from the unit.
        """
        if name not in self.partts:
            raise ValueError(f"unit has no '{name}' partition")
        self.cache['partitions'].pop(name)

    def getinterp(self, xnodes, ynodes):
        """Creates an interpolator on a mesh.

        Parameters
        ----------
        xnodes : flat-float-array
            x-coordinates of the interpolation nodes.
        ynodes : flat-float-array
            y-coordinates of the interpolation nodes.

        Returns
        -------
        TriInterp
            Callable interpolator.

        """
        return trinterp.getinterp(
            self.mesh, xnodes, ynodes
        )

    def massopr(self, is_lumped, add_constr):
        """Creates the mass operator from the base partition.

        Parameters
        ----------
        is_lumped : bool
            Creates a lumped mass operator, if True.
        add_constr : bool
            Adds constraints, if True.

        Returns
        -------
        MatrixFEM
            Mass operator as a matrix.

        """
        if is_lumped is True:
            return self.massopr_lumped(add_constr)
        return self.massopr_full(add_constr)

    def massopr_full(self, add_constr):
        return self.base.new_matrix(self.massmat, add_constr)

    def massopr_lumped(self, add_constr):
        matrix = self.base.new_matrix(self.massdiag, add_constr)
        return matrix.with_no_zeros()

    def constrproj(self):
        return femfactory.get_constr_proj(self)


class FEMUnitMaker:
    """Maker of FEM unit data.
    """

    def __init__(self):
        self.cache = {}

    def get_unit_data(self, mesh, anchors=None, mode='fem'):

        _ = self.make_mesh_aligned(mesh, anchors)
        _ = self.make_skeleton()
        _ = self.make_data(mode)

        return _

    def make_mesh_aligned(self, mesh, anchors):

        anchors = anchors or ()

        mesh0 = mesh
        mesh1 = mesh0.alignnodes(*anchors)

        mesh2 = mesh1.downvoids()
        mesh2 = mesh2.alignvoids()

        self.cache['meshes'] = {
            'root': mesh0,
            'alfa': mesh1,
            'beta': mesh2
        }

        return mesh2

    def make_skeleton(self):

        skel = skeleton.getskeleton(
            self.cache['meshes']['beta']
        )

        self.cache['skel'] = skel
        return skel

    def make_data(self, mode):

        mesh = self.make_unit_mesh()
        perm = self.make_data_perm()

        ij_v = self.make_ij_v_data(mode)

        meta = {
            'ij-stream': ij_v['ij-stream'],
            'v-streams': ij_v['v-streams'],
            'mesh2unit': perm,
            'unit-mode': mode
        }

        return {
            "mesh": mesh,
            "meta": meta
        }

    def make_unit_mesh(self):

        mesh = self.cache['meshes']['beta']
        mesh = mesh.twin()

        meta = {
            'voids': mesh.getvoids(),
            'loops': mesh.meshedge().getloops()
        }

        mesh = mesh.add_meta(meta)
        return mesh

    def make_ij_v_data(self, mode):

        skel = self.cache['skel']

        return {
            'ij-stream': ijstream.getstream(skel, mode),
            'v-streams': vstreams.getstreams(skel, mode)
        }

    def make_data_perm(self):

        mesh0 = self.cache['meshes']['root']
        mesh1 = self.cache['meshes']['alfa']
        mesh2 = self.cache['meshes']['beta']

        perm1 = mesh1.meta['nodes-permuter']
        perm2 = mesh2.meta['nodes-permuter']

        perm = perm1[perm2]

        meta = {
            'perm': perm,
            'perm-inv': np.argsort(perm)
        }

        return DataPermuter(mesh0, meta)


class DataPermuter:
    """Data permutator.

    Attributes
    ----------
    mesh : TriMesh
        Source triangle mesh.
    meta : dict
        Permutation metadata.

    """

    def __init__(self, mesh, meta):
        self.mesh = mesh
        self.meta = meta

    @property
    def perm(self):
        """Permutation of the sorce data.
        """
        return self.meta['perm']

    @property
    def perm_inv(self):
        """Inverse permutation.
        """
        return self.meta['perm-inv']

    def __call__(self, source_data):
        return self.permute(source_data)

    def permute(self, data):
        return np.copy(
            data[self.perm], order='C'
        )

    def permute_inv(self, data):
        return np.copy(
            data[self.perm_inv], order='C'
        )
