import numpy as np
import logging
from .sparsediffs import SparseDiffs, SparseValues


class Background:
    extension_sizes = [500, 5000]

    def __init__(self, linear_repr, path, config):
        self._fragment_length = config.fragment_length
        self._linear_repr = linear_repr
        self._path = path
        self._min_value = None

    def set_min_value(self, value):
        self._min_value = value

    def _get_pileup(self, positions, extension_size):
        start_ends = np.add.outer(np.array([-extension_size, extension_size]),
                                  positions)
        sparse_diffs = SparseDiffs.from_starts_and_ends(start_ends)
        sparse_diffs /= (extension_size*2/self._fragment_length)
        return sparse_diffs

    def create(self, reads):
        positions = [(read.node_ids[0], read.start) if read.direction > 0 else
                     (read.node_ids[-1], read.end-1) for read in reads]
        positions = np.array(positions)
        linear_positions = self._path.project_positions(positions)
        print(linear_positions)
        if self._min_value is None:
            self._min_value = linear_positions.size*self._fragment_length / self._path._distance_to_node[-1]
        logging.info("Using min value %s", self._min_value)
        max_pileup = None
        for extension_size in self.extension_sizes:
            pileup = self._get_pileup(positions, extension_size)
            if max_pileup is None:
                pileup.clip_min(self._min_value)
                max_pileup = pileup
            else:
                max_pileup = max_pileup.maximum(pileup)

        max_pileup._sanitize()
        pileup = max_pileup.get_sparse_values()
        nodes, offsets = self._path.back_project_positions(pileup.indices)
        linear_codes = self._linear_repr.get_linear_codes(nodes, offsets)
        return SparseValues(linear_codes, pileup.values)


class BackgroundFromSample(Background):
    extension_sizes = [5000]
