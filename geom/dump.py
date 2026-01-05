# -*- coding: utf-8 -*-
"""Makes geo files.
"""
import itertools as itr


def dump_path_map(pathmap, geopath, seeds):

    geofile = DumpPath.from_pathmap(pathmap).make_geofile(seeds)

    with open(geopath, encoding='utf-8', mode='w') as file:
        file.write(geofile)


class PathAgent:
    """Operator on a path-map.
    """

    def __init__(self, path=None):
        self.path = path
        self.cache = {}

    @classmethod
    def from_pathmap(cls, path):
        return cls(path)


class DumpPath(PathAgent):
    """Dumps a path-map.
    """

    def make_geofile(self, seeds: dict):

        points = self.make_points(seeds)
        lines = self.make_lines()
        surf = self.make_surf()

        return '\n\n'.join(
            [points, lines, surf]
        )

    def make_points(self, seeds):
        return DumpPoints.from_pathmap(self.path).make_points(seeds)

    def make_lines(self):
        return MakeLines().make_lines(self.path.size)

    def make_surf(self):
        return '\n'.join(
            [make_line_loop(self.path.size), make_planesurf()]
        )


class DumpPoints(PathAgent):
    """Dumps a path-map to points.
    """

    def make_points(self, seeds: dict):

        self.cache['seeds-dict'] = seeds

        coords = self.dump_coords()
        points = self.dump_points(coords)

        return '\n'.join(points)

    def dump_points(self, coords):
        return list(
            itr.starmap(self.dump_point, enumerate(coords))
        )

    def dump_point(self, i_d, coords):
        return f'Point({i_d + 1}) = ' + '{' + coords + '};'

    def dump_coords(self):
        return list(
            itr.starmap(self.dump_coord, self.coords_zipped)
        )

    def dump_coord(self, xpos, ypos, seed):
        return f'{xpos:.5f}, {ypos:.5f}, 0, {seed:.5f}'

    @property
    def seeds_dict(self):
        return self.cache['seeds-dict']

    @property
    def seeds_array(self):
        return [
            self.seeds_dict.get(color, 0) for color in self.path.colors
        ]

    @property
    def coords_zipped(self):

        x_pos = self.path.points.real
        y_pos = self.path.points.imag

        return zip(
            x_pos, y_pos, self.seeds_array
        )


class MakeLines:
    """Lines maker.
    """

    def make_lines(self, count):

        for_loop = self.make_for_loop(count)
        lastline = self.make_lastline(count)

        return '\n'.join([for_loop, lastline])

    def make_for_loop(self, count):

        head = self.make_for_head(count)
        body = self.make_for_body()
        tail = self.make_for_tail()

        return '\n'.join(
            [head, body, tail]
        )

    def make_for_head(self, count):
        return f'For i In {{1:{count-1}}}'

    def make_for_body(self):
        return '  Line(i) = {i, i+1};'

    def make_for_tail(self):
        return 'EndFor'

    def make_lastline(self, count):
        return f'Line({count}) = {{{count}, 1}};'


def make_line_loop(count):
    return f'Line Loop(1) = {{1:{count}}};'


def make_planesurf():
    return 'Plane Surface(1) = {1};'
