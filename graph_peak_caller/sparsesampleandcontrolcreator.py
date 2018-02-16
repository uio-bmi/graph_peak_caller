import numpy as np
from .sparsediffs import SparseDiffs
from .linearsnarls import LinearPileup
from .snarlmaps import LinearSnarlMap
import logging

class SparseControl:
    def __init__(self, linear_map, graph, extension_sizes, fragment_length, touched_nodes):
        self._linear_map = LinearSnarlMap.from_json_files(linear_map, graph)
        self._extension_sizes = extension_sizes
        self._fragment_length = fragment_length
        self._graph = graph
        self._min_value = None
        self._touched_nodes = touched_nodes

    def set_min_value(self, value):
        self._min_value = value

    def create(self, reads):
        mapped_reads = self._linear_map.map_interval_collection(reads)
        if self._min_value is None:
            self._min_value = mapped_reads.n_intervals*self._fragment_length / self._linear_map._length
        logging.info("Using min value %s", self._min_value)
        max_pileup = None
        # SparseDiffs([0], [self._min_value])
        for tmp_extension in self._extension_sizes:
            extension = tmp_extension // 2
            sparse_diffs = SparseDiffs.from_starts_and_ends(
                mapped_reads.extend_np(extension))
            sparse_diffs /= (extension*2/self._fragment_length)
            if max_pileup is None:
                sparse_diffs.clip_min(self._min_value)
                max_pileup = sparse_diffs
            else:
                max_pileup = max_pileup.maximum(sparse_diffs)
            # max_pileup = max_pileup.maximum(sparse_diffs)
        max_pileup._sanitize()
        lin_pileup = LinearPileup(
            max_pileup._indices,
            np.cumsum(max_pileup._diffs))
        lin_pileup.sanitize_indices()
        lin_pileup.sanitize_values()
        return lin_pileup.to_sparse_pileup(
            self._linear_map, self._touched_nodes, self._min_value)
