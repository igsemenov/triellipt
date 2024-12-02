# -*- coding: utf-8 -*-
"""Meshing tools.
"""
from triellipt.mesher import trigrids


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
        New mesh.

    Notes
    -----

    (i) Triangulation in terms of grid cell division:

    - "west-slope" — by west diagonals
    - "east-slope" — by east diagonals
    - "west-snake" — snake-like, starting from the west slope
    - "east-snake" — snake-like, starting from the east slope
    - "cross-wise" — by both diagonals

    """
    return trigrids.getgrid(
        int(xsize), int(ysize), slopes
    )
