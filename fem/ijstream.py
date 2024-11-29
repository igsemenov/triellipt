# -*- coding: utf-8 -*-
"""Stream of ij-coordinates.
"""
from triellipt.fem import ijstream_fvm
from triellipt.fem import ijstream_fem


STREAMERS = {
    'fvm': ijstream_fvm.getstream,
    'fem': ijstream_fem.getstream
}


def getstream(skel, mode):
    return STREAMERS[mode](skel)
