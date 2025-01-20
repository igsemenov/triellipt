# -*- coding: utf-8 -*-
"""Curves classes.
"""
from abc import ABC, abstractmethod
import numpy as np
from triellipt.geom import partt

CurveError = type(
    'CurveError', (Exception,), {}
)


class Curve(ABC):
    """Base class for parametric curves on [0, 1].
    """

    def getpath(self, *args):
        """Picks points on the curve.

        Parameters
        ----------
        args : *float
            Parameters in [0, 1].

        Returns
        -------
        flat-complex-array
            Points on the curve.

        """
        return self._getpath(np.array(args))

    @abstractmethod
    def _getpath(self, args):
        """Picks points on the curve for a given parameters array. 
        """

    def linspace(self, nparts):
        """Splits the curve uniformly in the parameter space.

        Parameters
        ----------
        nparts : int
            Number of intervals in the parameter space.

        Returns
        -------
        flat-complex-array
            The resulting polygonal path.

        """
        return self._getpath(
            np.linspace(0., 1., nparts + 1)
        )

    def partition(self, nparts, ratio=1):
        """Splits the curve into segments based on length.

        Parameters
        ----------
        nparts : int
            Number of segments in the partition.
        ratio : float = 1
            Ratio of segments lengths (last-to-first).

        Returns
        -------
        flat-complex-array
            The resulting polygonal chain.

        """

        ratio = ratio if ratio > 0. else 1.
        nseeds = partt.LenEstim(self).get_nseeds(rtol=1e-4, maxitr=10)

        if nseeds is None:
            raise CurveError(
                "cannot estimate the length, split the curve"
            )

        partter = partt.CurvePartt(self)
        return partter.get_partition(nparts, ratio, nseeds)

    def length(self, *, places=4, maxitr=10):
        """Estimates the curve length iteratively.

        Parameters
        ----------
        places : int = 4
            Number of decimal places to resolve.
        maxitr : int = 10
            Maximum number of iterations.

        Returns
        -------
        float
            Length estimate.

        """

        rtol = pow(10, -places)
        length = partt.LenEstim(self).get_length(rtol, maxitr)

        if length is not None:
            return length

        raise CurveError(
            "cannot estimate the length, maximum iteration reached"
        )

    def is_curve(self):
        return isinstance(
            self, (Line, EllipArc, Bezier2, Bezier3)
        )


class Line(Curve):
    """Line between two points.
    """

    def __init__(self, startpoint, endpoint):
        self._startpoint: complex = startpoint
        self._endpoint: complex = endpoint

    @property
    def startpoint(self):
        return self._startpoint

    @property
    def endpoint(self):
        return self._endpoint

    def _getpath(self, args):
        return (1. - args) * self.startpoint + args * self.endpoint


class EllipArc(Curve):
    """Elliptic arc.
    """

    def __init__(self, center, axes, phis, tilt=0):

        self.center: complex = center

        self.axes: tuple = axes  # major-minor-axes
        self.phis: tuple = phis  # start-end-arc-angles
        self.tilt: float = tilt  # arc-rotation-angle

    def _getpath(self, args):

        phi = self.linphi(args)

        axsum = self.axes[0] + self.axes[1]
        axdif = self.axes[0] - self.axes[1]

        phase1 = 0.5 * axsum * np.exp(+1j * phi)
        phase2 = 0.5 * axdif * np.exp(-1j * phi)

        path = self.center + phase1 + phase2

        if not self.tilt:
            return path
        return np.exp(1j * self.tilt) * path

    def linphi(self, args):
        return (1. - args) * self.phis[0] + args * self.phis[1]

    @property
    def startpoint(self):
        return self.getpath(0.)

    @property
    def endpoint(self):
        return self.getpath(1.)


class Bezier(Curve):
    """ABC for Bezier curves.
    """

    def __init__(self, *control_points):

        self.points: tuple[complex] = control_points

        self.westchord = self.getwestchord()  # west-bezier-base
        self.eastchord = self.geteastchord()  # east-bezier-base

    def _getpath(self, args):

        westchord = self.westchord.getpath(*args)
        eastchord = self.eastchord.getpath(*args)

        return westchord * (1. - args) + eastchord * args

    @abstractmethod
    def getwestchord(self):
        """Returns the west curve in a linear combination.
        """

    @abstractmethod
    def geteastchord(self):
        """Returns the east curve in a linear combination.
        """

    @property
    def startpoint(self):
        return self.points[0]

    @property
    def endpoint(self):
        return self.points[-1]


class Bezier2(Bezier):
    """Quadratic Bezier curve.
    """

    def __init__(self, point0, point1, point2):
        super().__init__(point0, point1, point2)

    def getwestchord(self):
        return Line(
            self.points[0], self.points[1]
        )

    def geteastchord(self):
        return Line(
            self.points[1], self.points[2]
        )


class Bezier3(Bezier):
    """Cubic Bezier curve.
    """

    def __init__(self, point0, point1, point2, point3):
        super().__init__(point0, point1, point2, point3)

    def getwestchord(self):
        return Bezier2(
            self.points[0], self.points[1], self.points[2]
        )

    def geteastchord(self):
        return Bezier2(
            self.points[1], self.points[2], self.points[3]
        )
