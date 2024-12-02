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
