# -*- coding: utf-8 -*-
"""Triangle fronts.
"""
import numpy as np
from triellipt.fem import skeleton
from triellipt.utils import tables


class TriFront:
    """Front of triangles.

    Properties
    ----------

     Name       | Description
    ------------|----------------------------------------
    `trinums`   | Indices of the front-facing triangles.
    `voidnums`  | Indices of void triangles in the front.

    """

    def __init__(self, unit=None, data=None):
        self.unit = unit
        self.data = data

    def update_data(self, new_data):
        return self.__class__(
            self.unit, new_data.copy('C')
        )

    def subfront(self, *inds):
        return self.update_data(self.data[:, inds])

    @property
    def size(self):
        return self.data.shape[1]

    @property
    def trinums(self):
        return self.data[0, :]

    @property
    def voidnums(self):
        return self.data[1, :]

    @property
    def trinums_voids(self):
        return self.data[1, :]

    @property
    def centrs2d(self):
        return _unpack_complex(self.centrs_complex)

    @property
    def centrs2d_voids(self):
        return _unpack_complex(
            self.centrs_complex_voids
        )

    @property
    def centrs_complex(self):
        return self.unit.mesh.centrs_complex[self.trinums]

    @property
    def centrs_complex_voids(self):
        return self.unit.mesh.centrs_complex[self.trinums_voids]

    @property
    def ranks(self):
        return self.getranks()

    def getranks(self):

        _, counts = np.unique(
            self.trinums, return_counts=True
        )

        return np.unique(counts).tolist()

    def atranks(self):

        inds = dict(
            _spec_col1d(self.data[[0], :].T)
        )

        data = {
            r: self.subfront(*i) for r, i in inds.items()
        }

        return data

    def atrank(self, rank):
        """Selects the front with the specified rank.
        """
        return self.atranks()[rank]

    def angles(self):
        """Computes the orientation angles of the front.

        Returns
        -------
        flat-float-array
            Angles between the voids and the front centroids.

        """
        return np.angle(
            self.centrs_complex - self.centrs_complex_voids
        )

    def scales(self):
        """Computes the normalized front scales.

        Returns
        -------
        flat-float-array
            Normalized distances between the front centroids and voids.

        """

        scales = np.abs(
            self.centrs_complex - self.centrs_complex_voids
        )

        return _normalize(scales)

    def filter_by_mask(self, mask):
        """Filters the front by the mask.

        Parameters
        ----------
        mask : function
            Boolean mask `(x, y) ` on the triangles centroids.

        Returns
        -------
        TriFront
            New front.

        """
        return self.update_data(
            self.data[:, mask(*self.centrs2d)]
        )

    def filter_by_angle(self, angmin, angmax):
        """Filters the front by the orientation angles.
        """

        angs = self.angles()

        mask = np.logical_and(
            angs >= angmin, angs <= angmax
        )

        return self.update_data(self.data[:, mask])

    def filter_by_scale(self, minval, maxval):
        """Filters the front by the scales.
        """

        scales = self.scales()

        mask = np.logical_and(
            scales >= minval, scales <= maxval
        )

        return self.update_data(self.data[:, mask])


class TriFrontFine(TriFront):
    """Front of fine triangles.
    """

    @classmethod
    def from_unit(cls, unit):
        return MakerFrontFine.from_unit(unit).get_front()

    @property
    def trinums_wests(self):
        return self.data[2, :]

    @property
    def trinums_easts(self):
        return self.data[3, :]

    def suptriu(self):
        """Returns the associated supertriangulation.
        """

        if self.size == 0:
            return None

        _ = self.unit.mesh.supertriu()

        return _.atcores(
            *np.unique(self.trinums)
        )


class TriFrontCoarse(TriFront):
    """Front of coarse triangles.
    """

    @classmethod
    def from_unit(cls, unit):
        return MakerFrontCoarse.from_unit(unit).get_front()


class UnitAgent:
    """Operator on an AMR unit.
    """

    def __init__(self, unit):
        self.unit = unit
        self.meta = self.fetch_meta()
        self.cache = {}

    @classmethod
    def from_unit(cls, unit):
        return cls(unit)

    def fetch_meta(self):
        return {}


class MakerFrontCoarse(UnitAgent):
    """Finds a front of coarse triangles.
    """

    def fetch_meta(self):

        suptri = self.unit.mesh.supertriu()

        return {
            'supvoids': suptri.supvoids()
        }

    @property
    def supvoids(self):
        return self.meta['supvoids']

    def get_front(self):

        dat = self.make_front_data()
        obj = self.make_front_object(dat)

        return obj

    def make_front_data(self):

        if self.supvoids is None:
            return np.array([[]]*2).astype(int)

        data = _pack_rows(
            self.supvoids.trinums1,
            self.supvoids.trinums
        )

        return data.astype(int)

    def make_front_object(self, data):
        return TriFrontCoarse(self.unit, data)


class MakerFrontFine(UnitAgent):
    """Finds a front of fine triangles.
    """

    def fetch_meta(self):

        skel = skeleton.getskeleton(self.unit.mesh)

        return {
            'voidsmap': skel.voidsmap
        }

    @property
    def voidsmap(self):
        return self.meta['voidsmap']

    def get_front(self):

        dat = self.make_front_data()
        obj = self.make_front_object(dat)

        return obj

    def make_front_data(self):

        if self.voidsmap is None:
            return np.array([[]]*4).astype(int)

        cores = self.voidsmap.trinums[2::4]
        voids = self.voidsmap.trinums[0::4]
        wests = self.voidsmap.trinums[1::4]
        easts = self.voidsmap.trinums[3::4]

        return _pack_rows(
            cores, voids, wests, easts
        )

    def make_front_object(self, data):
        return TriFrontFine(self.unit, data)


def _pack_rows(*rows):
    return np.vstack(rows).copy('C')


def _unpack_complex(argz):
    return np.vstack(
        [argz.real, argz.imag]
    )


def _normalize(data):
    return data / np.amax(data)


def _spec_col1d(col1d):

    colmap = tables.maptable(col1d)
    rankset = np.unique(colmap.vals_ranks)

    for rank in rankset:
        _, rows, _ = colmap.atrank(rank)
        yield (rank, rows)
