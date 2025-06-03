# -*- coding: utf-8 -*-
"""Meshing tools.
"""
from triellipt.mesher import trigrids
from triellipt.mesher import trilattices


def trigrid(xsize, ysize, slopes):
    """Creates a triangle grid.

    Parameters
    ----------
    xsize : int
        Number of nodes in x-direction.
    ysize : int
        Number of nodes in y-direction.
    slopes : str
        Controls orientation of triangles (i).

    Returns
    -------
    TriMesh
        Resulting triangle mesh.

    Notes
    -----

    (i) Triangulation in terms of grid cell division:

    Key           | Description
    --------------|----------------------------------------------
    "west-slope"  | By west diagonals
    "east-slope"  | By east diagonals
    "cross-wise"  | By both diagonals
    "west-snake"  | Snake-like, starting from the west slope
    "east-snake"  | Snake-like, starting from the east slope

    """
    return trigrids.getgrid(
        int(xsize), int(ysize), slopes
    )


def trilattice(xsize, ysize, close=False):
    """Creates a lattice of equilateral triangles.

    Parameters
    ----------
    xsize : int
        Number of nodes in x-direction.
    ysize : int
        Number of nodes in y-direction.
    close : bool = False
        Closes the lattice sides, if True.

    Returns
    -------
    TriMesh
        Resulting triangle mesh.

    """
    return trilattices.get_lattice(xsize, ysize, close)
