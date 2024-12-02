# -*- coding: utf-8 -*-
"""Public functions.
"""
import math
from triellipt.geom import curves as curves_
from triellipt.geom import loop
from triellipt.geom import cycle


def line(startpoint, endpoint):
    """Creates a line between two points.

    Parameters
    ----------
    startpoint : complex
        Start point.
    endpoint : complex
        End point.

    Returns
    -------
    Line
        Line as a curve-like object.

    """
    return curves_.Line(startpoint, endpoint)


def elliparc(center, axes, phis, tilt=0):
    """Creates an elliptic arc.

    Parameters
    ----------
    center : complex
        Ellipse center.
    axes : (float, float)
        Major and minor ellipse axes.
    phis : (float, float)
        Start and end arc angles.
    tilt : float = 0
        Ellipse rotation angle around the center. 

    Returns
    -------
    EllipArc
        Elliptic arc as a curve-like object.

    """
    return curves_.EllipArc(center, axes, phis, tilt)


def bezier2(point0, point1, point2):
    """Creates a quadratic Bezier curve.

    Parameters
    ----------
    point0 : complex
        1st control point (startpoint).
    point1 : complex
        2nd control point.
    point2 : complex
        3rd control point (endpoint).

    Returns
    -------
    Bezier2
        Curve-like object.

    """
    return curves_.Bezier2(point0, point1, point2)


def bezier3(point0, point1, point2, point3):
    """Creates a cubic Bezier curve.

    Parameters
    ----------
    point0 : complex
        1st control point (startpoint).
    point1 : complex
        2nd control point.
    point2 : complex
        3rd control point.
    point3 : complex
        4th control point (endpoint).

    Returns
    -------
    Bezier3
        Curve-like object.

    """
    return curves_.Bezier3(
        point0, point1, point2, point3
    )


def makeloop(*curves):
    """Creates a loop of connected curves.

    Parameters
    ----------
    curves : *one-of-the-curves
        Sequence of connected curves.

    Returns
    -------
    CurvesLoop
        Loop of connected curves.

    """
    return loop.CurvesLoop.from_curves(*curves)


def makerect(corner, dims):
    """Creates a rectangle as a loop.

    Parameters
    ----------
    corner : complex
        South-west rectangle corner.
    dims : (float, float)
        Width and height of the rectangle.

    Returns
    -------
    CurvesLoop
        Loop of rectangle sides (south-east-north-west).

    """

    width, height = dims

    dispx = width
    dispy = height * 1j

    points = [
        corner,
        corner + dispx,
        corner + dispx + dispy,
        corner + dispy
    ]

    lines = [
        line(points[0], points[1]),
        line(points[1], points[2]),
        line(points[2], points[3]),
        line(points[3], points[0])
    ]

    return makeloop(*lines)


def makeellip(center, axes, tilt=0):
    """Creates a closed ellipse as a loop.

    Parameters
    ----------
    center : complex
        Ellipse center.
    axes : (float, float)
        Major and minor axes.
    tilt : float = 0
        Ellipse rotation angle.

    Returns
    -------
    CurvesLoop
        Ellipse as a single-curve loop.

    """

    curve = elliparc(
        center, axes, (0., 2.*math.pi), tilt
    )

    return makeloop(curve)


def makecycle(path):
    """Creates a cyclic path.

    Parameters
    ----------
    path : flat-complex-array
        Input path.

    Returns
    -------
    CycPath
        Path closed to a cycle.

    """
    return cycle.CycPath.from_path(path)
