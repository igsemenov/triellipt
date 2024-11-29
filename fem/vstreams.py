# -*- coding: utf-8 -*-
"""Stream of matrix entries.
"""
from triellipt.fem import vstreams_fvm
from triellipt.fem import vstreams_fem

STREAMERS = {
    'fvm': vstreams_fvm.getstreams,
    'fem': vstreams_fem.getstreams
}


def getstreams(skel, mode):
    return STREAMERS[mode](skel)
