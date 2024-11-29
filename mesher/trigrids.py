# -*- coding: utf-8 -*-
"""Triangle grids.
"""
from abc import ABC, abstractmethod
import numpy as np
from triellipt.trimesh import trimesh_

TriGridsError = type(
    'TriGridsError', (Exception,), {}
)


def getgrid(xsize, ysize, key):

    xsize = max(xsize, 2)
    ysize = max(ysize, 2)

    if key not in TRIGRIDS:
        raise TriGridsError(
            f"got undefined grid type '{key}'"
        )

    return TRIGRIDS[key](xsize, ysize).get_trimesh()


class TriGridder(ABC):
    """Triangle grid maker.
    """

    def __init__(self, xsize, ysize):
        self.xsize = xsize
        self.ysize = ysize

    def get_trimesh(self):
        points, triangs = self.gen_meshdata()
        return trimesh_.TriMesh.from_data(points, triangs)

    def gen_meshdata(self):
        yield self.get_points()
        yield self.get_triangs()

    @abstractmethod
    def get_points(self):
        """Creates the mesh points.

        Returns
        -------
        complex-flat-array
            Mesh points.

        """

    def get_triangs(self):
        """Creates the triangles table.

        Returns
        -------
        3-column-int-table
            Table of mesh triangles.

        """

        cells = self.make_cells()
        triangs = self.from_cells(cells)

        return triangs

    @abstractmethod
    def make_cells(self):
        """Returns table of grid cells.
        """

    @abstractmethod
    def from_cells(self, cells):
        """Makes triangles from grid cells.
        """

    def pack_triangs(self, *tables):

        tris_flat = np.hstack(tables).flatten()

        triangs = np.reshape(
            tris_flat, (tris_flat.size // 3, 3)
        )

        return triangs.copy('C')


class TriOneSlope(TriGridder):
    """Grid via one-slope splitting of cells.
    """

    TRIS_EAST = (
        (0, 1, 2),
        (2, 3, 0)
    )

    TRIS_WEST = (
        (1, 2, 3),
        (3, 0, 1)
    )

    def get_points(self):
        return self.numgrid.get_points()

    def make_cells(self):
        """Creates 4-column table of cells corners.
        """
        return self.numgrid.get_cells_m4()

    def gen_snake_masks(self):

        mask = self.numgrid.snake_mask_flat()

        yield mask
        yield np.logical_not(mask)

    @property
    def numgrid(self):
        return NumGrid(self.xsize, self.ysize)


class TriEastSlope(TriOneSlope):
    """Grid via splitting of cells by the east diagonal.

    Triangles
    ---------
      0-1-2
      2-3-0

    """

    def from_cells(self, cells):

        trisone = cells[:, self.TRIS_EAST[0]]
        tristwo = cells[:, self.TRIS_EAST[1]]

        return self.pack_triangs(trisone, tristwo)


class TriWestSlope(TriOneSlope):
    """Grid via splitting of cells by the west diagonal.

    Triangles
    ---------
      1-2-3
      3-0-1

    """

    def from_cells(self, cells):

        trisone = cells[:, self.TRIS_WEST[0]]
        tristwo = cells[:, self.TRIS_WEST[1]]

        return self.pack_triangs(trisone, tristwo)


class TriSnake(TriOneSlope):
    """Grid via splitting by diagonals in a snake mode.
    """

    def from_cells(self, cells):

        mask1, mask2 = self.gen_snake_masks()
        path11, path12, path21, path22 = self.gen_paths()

        tris11 = cells[mask1, :]
        tris12 = cells[mask1, :]

        tris21 = cells[mask2, :]
        tris22 = cells[mask2, :]

        tris11 = tris11[:, path11]
        tris12 = tris12[:, path12]

        tris21 = tris21[:, path21]
        tris22 = tris22[:, path22]

        trisone = np.vstack([tris11, tris21])
        tristwo = np.vstack([tris12, tris22])

        return self.pack_triangs(trisone, tristwo)

    def gen_paths(self):

        paths1 = self.get_paths_one()
        paths2 = self.get_paths_two()

        yield paths1[0]
        yield paths1[1]
        yield paths2[0]
        yield paths2[1]

    @abstractmethod
    def get_paths_one(self):
        """Returns paths for the east-part.
        """

    @abstractmethod
    def get_paths_two(self):
        """Returns paths for the west-part.
        """


class TriEastSnake(TriSnake):
    """Grid via splitting of cells in the east-snake mode.
    """

    def get_paths_one(self):
        return self.TRIS_EAST

    def get_paths_two(self):
        return self.TRIS_WEST


class TriWestSnake(TriSnake):
    """Grid via splitting of cells in the west-snake mode.
    """

    def get_paths_one(self):
        return self.TRIS_WEST

    def get_paths_two(self):
        return self.TRIS_EAST


class TriCrossWise(TriGridder):
    """Grid via crosswise splitting of cells.

    Triangles
    ---------
    4-3-0
    4-0-1
    4-1-2
    4-2-3

    """

    def get_points(self):

        points1 = self.numgrid_base.get_points()
        points2 = self.numgrid_core.get_points()

        points = [
            points1, points2 + (0.5 + 0.5 * 1j)
        ]

        return np.hstack(points)

    def make_cells(self):
        return self.numgrid_base.get_cells_m5()

    def from_cells(self, cells):

        tris1 = cells[:, [4, 3, 0]]
        tris2 = cells[:, [4, 0, 1]]
        tris3 = cells[:, [4, 1, 2]]
        tris4 = cells[:, [4, 2, 3]]

        return self.pack_triangs(
            tris1, tris2, tris3, tris4
        )

    @property
    def numgrid_base(self):
        return NumGrid(self.xsize, self.ysize)

    @property
    def numgrid_core(self):
        return self.numgrid_base.coarsen()


class NumGrid:
    """2D grid of indices.

    Grid cell
    ---------
      3 — 2
      | 4 |
      0 — 1

    """

    def __init__(self, xsize, ysize):

        self.grid = np.reshape(
            np.arange(xsize*ysize), (xsize, ysize)
        )

        self.xsize = xsize
        self.ysize = ysize

    @property
    def offset(self):
        """Maximum number in the grid.
        """
        return self.xsize * self.ysize

    @property
    def subgrid(self):
        return self.coarsen()

    def coarsen(self):
        return self.__class__(
            self.xsize - 1, self.ysize - 1
        )

    def get_points(self):

        xpos, ypos = map(
            np.ravel, np.mgrid[:self.xsize, :self.ysize]
        )

        return xpos + 1j * ypos

    def get_cells_m4(self):
        return _packgen(
            self.gen_cells_m4()
        )

    def get_cells_m5(self):
        return _packgen(
            self.gen_cells_m5()
        )

    def gen_cells_m4(self):
        yield self.grid[:-1, :-1].flatten()
        yield self.grid[1::, :-1].flatten()
        yield self.grid[1::, 1::].flatten()
        yield self.grid[:-1, 1::].flatten()

    def gen_cells_m5(self):
        for data in self.gen_cells_m4():
            yield data
        yield self.get_cores()

    def get_cores(self):
        subgrid = self.subgrid.grid
        return subgrid.flatten() + self.offset

    def snake_mask_flat(self):
        return _snake_mask(
            self.xsize - 1, self.ysize - 1
        )


def _packgen(genobj):
    """Packs generator of arrays to a table.
    """
    return np.vstack(list(genobj)).T.copy('C')


def _snake_mask(isize, jsize):

    mask = np.full(
        (isize, jsize), False
    )

    mask[0::2, 0::2] = True
    mask[1::2, 1::2] = True

    return mask.flatten()


TRIGRIDS = {
    'west-slope': TriWestSlope,
    'east-slope': TriEastSlope,
    'west-snake': TriWestSnake,
    'east-snake': TriEastSnake,
    'cross-wise': TriCrossWise
}
