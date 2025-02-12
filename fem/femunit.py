# -*- coding: utf-8 -*-
"""Public FEM unit.
"""
import numpy as np
from triellipt.fem import (
    skeleton,
    ijstream,
    vstream,
    femfactory,
    femvector,
    femoprs,
    fempartt,
    trinterp
)


def getunit(mesh, anchors=None):
    """Creates a FEM computing unit.

    Parameters
    ----------
    mesh : TriMesh
        Input triangle mesh.
    anchors : tuple = None
        Nodes numbers to synchronize the mesh boundary.

    Returns
    -------
    FEMUnit
        FEM computing unit.

    """
    return FEMUnit.from_mesh(mesh, anchors)


class FEMData:
    """Base FEM data.
    """

    def __init__(self, mesh=None, meta=None):
        self.mesh = mesh
        self.meta = meta
        self.cache = {}

    @classmethod
    def from_mesh(cls, mesh, anchors=None):
        return FEMUnitMaker().get_unit(mesh, anchors)

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
        return self.femoprs['diff_1x'].data

    @property
    def diff_1y(self):
        return self.femoprs['diff_1y'].data

    @property
    def diff_2x(self):
        return self.femoprs['diff_2x'].data

    @property
    def diff_2y(self):
        return self.femoprs['diff_2y'].data

    @property
    def femoprs(self):
        return self.meta['femmatrix']['v-stream']

    @property
    def ij_stream(self):
        return self.meta['femmatrix']['ij-stream']

    @property
    def ij_tuple(self):
        return (
            self.ij_stream.rownums, self.ij_stream.colnums
        )

    @property
    def grad(self):
        """Gradient operator.
        """

        if 'grad' in self.cache:
            return self.cache['grad']

        self.cache['grad'] = femoprs.getgrad(self.mesh)
        return self.grad

    @property
    def perm(self):
        """Nodes permutation in the mesh-to-unit transition.
        """
        return self.meta['data-perm']


class FEMUnit(FEMRoot):
    """FEM computing unit.

    Properties
    ----------

    FEM operators as data-streams:

    Name        | Description
    ------------|---------------------
    `massmat`   | Mass-matrix
    `massdiag`  | Mass-matrix lumped
    `diff_1y`   | 1st-y derivative
    `diff_1x`   | 1st-x derivative
    `diff_2y`   | 2nd-y derivative
    `diff_2x`   | 2nd-x derivative

    FEM operators as FEM matrices:

    Name           | Description
    -------------- |----------------------------------------
    `laplace_fem`  | Laplace operator
    `massmat_fem`  | Mass-matrix (no constraints)
    `massmat_amr`  | Mass-matrix (with constraints)
    `massdiag_fem` | Mass-matrix lumped (no constraints)
    `massdiag_amr` | Mass-matrix lumped (with constraints)

    Others:

    Name      | Description
    ----------|----------------------------------
    `grad`    | Gradient operator.
    `perm`    | Mesh-to-unit nodes permutation.
    `loops`   | List of the mesh loops.

    """

    def fem_factory(self, add_constraints=True):
        """Creates a factory of FEM matrices.

        Parameters
        ----------
        add_constraints : bool = True
            If True, forces constraints to be included in the matrix.

        Returns
        -------
        FEMFactory
            Resulting factory of FEM matrices.

        """
        return femfactory.FEMFactory.from_unit(self, add_constraints)

    def new_vector(self):
        """Returns a new FEM vector.
        """
        return femvector.VectorFEM.from_unit(self)

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
            Interpolator object.

        """
        return trinterp.getinterp(self, xnodes, ynodes)

    def partition(self, loopnum, splitang, parttspec=None):
        """Creates the FEM unit partition.

        Parameters
        ----------
        loopnum : int
            Number of the mesh loop to partition.
        splitang : float
            Threshold angle for the loop splitting.
        parttspec : list = None
            List of `(edge1, edge2, operation)` triplets (a).

        Returns
        -------
        Partition
            Partition object.

        Notes
        -----

        (a) Color-operations on the edges:

        - "l" — left-repainting
        - "r" — right-repainting
        - "s" — switching-of-colors

        """
        return fempartt.getpartt(
            self, loopnum, splitang, parttspec
        )

    @property
    def laplace_fem(self):
        return self.fct1.feed_data(
            self.diff_2x + self.diff_2y
        )

    @property
    def massmat_fem(self):
        return self.fct0.feed_data(self.massmat)

    @property
    def massdiag_fem(self):
        return self.fct0.feed_data(self.massdiag).with_no_zeros()

    @property
    def massmat_amr(self):
        return self.fct1.feed_data(self.massmat)

    @property
    def massdiag_amr(self):
        return self.fct1.feed_data(self.massdiag).with_no_zeros()

    @property
    def constrproj(self):
        return femfactory.get_constr_proj(self)

    @property
    def fct0(self):
        """Factory with no constraints.
        """

        if 'fct0' in self.cache:
            return self.cache['fct0']

        self.cache['fct0'] = self.fem_factory(0)
        return self.fct0

    @property
    def fct1(self):
        """Factory with constraints.
        """

        if 'fct1' in self.cache:
            return self.cache['fct1']

        self.cache['fct1'] = self.fem_factory(1)
        return self.fct1


class FEMUnitMaker:
    """Maker of FEM unit data.
    """

    def __init__(self):
        self.cache = {}

    def get_unit(self, mesh, anchors=None):

        _ = self.make_mesh_aligned(mesh, anchors)
        _ = self.make_skeleton()
        _ = self.make_unit()

        return _

    def make_mesh_aligned(self, mesh, anchors):

        anchors = anchors or ()

        mesh0 = mesh
        mesh1 = mesh0.alignnodes(*anchors)
        mesh2 = mesh1.alignvoids()

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

    def make_unit(self):

        mesh = self.make_unit_mesh()
        ij_v = self.make_ij_v_data()
        perm = self.make_data_perm()

        meta = {
            'femmatrix': ij_v,
            'data-perm': perm
        }

        return FEMUnit(mesh, meta)

    def make_unit_mesh(self):

        mesh = self.cache['meshes']['beta']
        mesh = mesh.twin()

        meta = {
            'voids': mesh.getvoids(),
            'loops': mesh.meshedge().getloops()
        }

        mesh = mesh.add_meta(meta)
        return mesh

    def make_ij_v_data(self):

        skel = self.cache['skel']

        return {
            'ij-stream': ijstream.getstream(skel),
            'v-stream': vstream.getstream(skel)
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
            'perm-inv':  np.argsort(perm)
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
