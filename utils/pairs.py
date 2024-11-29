# -*- coding: utf-8 -*-
"""Pairing of integers.
"""
import itertools as itr
import numpy as np


def sympaired(ipos, jpos):
    """Symmetric pairing of two integers.
    """

    ipos = ipos.astype(np.int64)
    jpos = jpos.astype(np.int64)

    upair = ipos + jpos * jpos
    lpair = jpos + ipos * ipos

    return np.where(
        ipos < jpos, upair, lpair
    )


def szupaired(ipos, jpos):
    """Szudzik pairing of two integers.
    """

    ipos = ipos.astype(np.int64)
    jpos = jpos.astype(np.int64)

    upair = ipos + jpos * jpos
    lpair = jpos + ipos * ipos + ipos

    return np.where(
        ipos < jpos, upair, lpair
    )


def szuunpaired(argz):
    """Szudzik unpairing function.
    """

    sqrtz = np.array(
        np.floor(np.sqrt(argz)), dtype=np.int64
    )

    dist1 = argz - sqrtz * sqrtz
    dist2 = argz - sqrtz * sqrtz - sqrtz

    argx = np.where(dist2 < 0, dist1, sqrtz)
    argy = np.where(dist2 < 0, sqrtz, dist2)

    return argx, argy


def paircols(table):
    """Pairing of table columns.

    Parameters
    ----------
    table : n-column-int-table
        Table to process.

    Returns
    -------
    n-column-int64-table
        Table of paired columns.

    """

    colpairs = [
        sympaired(i, j) for i, j in itr.pairwise(table.T)
    ]

    return np.vstack(colpairs).T.copy('C')
