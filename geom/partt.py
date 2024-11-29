# -*- coding: utf-8 -*-
"""Curve partition.
"""
import numpy as np


class CurveAgent:
    """Operator on a curve.
    """

    def __init__(self, curve):
        self.curve = curve

    def lenfunc_from_linspace(self, nseeds):

        path = self.curve.linspace(nseeds)
        pins = np.abs(np.diff(path))

        return np.add.accumulate(pins)

    def length_from_linspace(self, nseeds):

        path = self.curve.linspace(nseeds)
        pins = np.abs(np.diff(path))

        return np.sum(pins)


class CurvePartt(CurveAgent):
    """Curve partitioner.
    """

    def get_partition(self, nparts, ratio, nseeds):
        """Returns the curve partition as a complex-flat-array.
        """

        args = self.find_args(nparts, ratio, nseeds)
        path = self.from_args(args)

        return path

    def from_args(self, args):
        return self.curve.getpath(*args)

    def find_args(self, nparts, ratio, nseeds):

        data = self.create_lenmap(nparts, ratio, nseeds)
        args = self.invert_lenmap(data)

        return args

    def create_lenmap(self, nparts, ratio, nseeds):

        func_we_know = self.lenfunc_we_know(nseeds)
        func_we_want = self.lenfunc_we_want(nparts, ratio, nseeds)

        return {
            'func-we-know': func_we_know,
            'func-we-want': func_we_want
        }

    def lenfunc_we_know(self, nseeds):
        return self.lenfunc_from_linspace(nseeds)

    def lenfunc_we_want(self, nparts, ratio, nseeds):

        usteps = self.usteps_from_params(nparts, ratio)
        length = self.length_from_linspace(nseeds)

        return np.add.accumulate(
            usteps * (length / np.sum(usteps))
        )

    def usteps_from_params(self, nparts, ratio):
        return 1. + np.linspace(0., 1., nparts) * (ratio - 1.)

    def invert_lenmap(self, lenmap):

        func_we_know = lenmap['func-we-know']
        func_we_want = lenmap['func-we-want']

        func_we_know = _normalize_lenfunc(func_we_know)
        func_we_want = _normalize_lenfunc(func_we_want)

        args_we_know = np.linspace(
            0., 1., func_we_know.size
        )

        args_we_want = np.interp(
            func_we_want, func_we_know, args_we_know
        )

        return np.r_[
            0., args_we_want[1:-1], 1.
        ]


def _normalize_lenfunc(data):
    return np.r_[
        0., data / data[-1]
    ]


class LenEstim(CurveAgent):
    """Curve length estimator.
    """

    def get_length(self, rtol, maxitr) -> float:
        """Returns the length estimate.
        """

        length, _ = self.estimate(
            getlen=self.length_from_linspace, rtol=rtol, maxitr=maxitr
        )

        return length

    def get_nseeds(self, rtol, maxitr) -> int:
        """Number of seeds we need to get the length estimate.
        """

        _, nseeds = self.estimate(
            getlen=self.length_from_linspace, rtol=rtol, maxitr=maxitr
        )

        return nseeds

    def estimate(self, getlen, rtol, nseeds=5, maxitr=10):

        prev = getlen(nseeds)

        for _ in range(maxitr):

            nseeds *= 2

            curr = getlen(nseeds)
            rerr = abs(curr/prev - 1.)

            if rerr < rtol:
                return curr, nseeds
            prev = curr

        return None, None
