# -*- coding: utf-8 -*-
"""Refiner dumper.
"""
import numpy as np
from scipy import sparse as sp


class DumpRefiner:
    """Refiner dumper.
    """

    def __init__(self, refiner):
        self.ref = refiner
        self.cache = {}

    @property
    def master_size(self):
        return self.cache['master-size']

    @property
    def slaves_size(self):
        return self.cache['slaves-size']

    def get_matrix(self):
        """Dumps the refiner to a matrix.
        """

        i_stream = self.make_i_stream()
        j_stream = self.make_j_stream()
        v_stream = self.make_v_stream()

        coo = sp.coo_array(
            (v_stream, (i_stream, j_stream))
        )

        return coo.tocsr()

    def make_v_stream(self):

        v_master = np.full(self.master_size, 1.0)
        v_slaves = np.full(self.slaves_size, 0.5)

        return np.hstack(
            [v_master, v_slaves]
        )

    def make_i_stream(self):

        i_master = self.ref.master_nodes

        i_slaves = np.repeat(
            self.ref.slaves_nodes, 2
        )

        self.cache = {
            'master-size': i_master.size,
            'slaves-size': i_slaves.size
        }

        return np.hstack(
            [i_master, i_slaves]
        )

    def make_j_stream(self):

        j_master = self.ref.nodes_images[self.ref.master_nodes, 0]
        j_slaves = self.ref.nodes_images[self.ref.slaves_nodes, :]

        return np.hstack(
            [j_master, j_slaves.flatten()]
        )
