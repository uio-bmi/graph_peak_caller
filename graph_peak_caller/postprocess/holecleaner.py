import logging
import numpy as np

from .segmentanalyzer import SegmentAnalyzer
from .graphs import DividedLinegraph, DummyTouched
from ..sparsediffs import SparseValues


class HolesCleaner:
    def __init__(self, graph, sparse_values, max_size, touched_nodes=None):
        sparse_values.indices = sparse_values.indices.astype("int")
        self._graph = graph
        self._node_indexes = graph.node_indexes
        self._sparse_values = sparse_values
        self._holes = self.get_holes()
        self._node_ids = self.get_node_ids()
        self._max_size = max_size
        self._kept = []
        self._touched_nodes = touched_nodes
        if touched_nodes is None:
            self._touched_nodes = DummyTouched()

    def get_holes(self):
        start_idx = 0
        if self._sparse_values.values[0] != 0:
            start_idx += 1
        end_idx = self._sparse_values.indices.size
        self._last_node = self._node_indexes.size
        if self._sparse_values.values[-1] == 0:
            self._last_node = np.searchsorted(
                self._node_indexes, self._sparse_values.indices[-1])
            end_idx -= 1
        return self._sparse_values.indices[start_idx:end_idx].reshape(
            (end_idx-start_idx)//2, 2)

    def build_sparse_values(self, holes):
        indices = np.sort(holes.ravel())
        values = np.arange(indices.size) % 2
        if (not indices.size) or indices[0] != 0:
            indices = np.r_[0, indices]
            values = np.r_[1, values]
        if self._sparse_values.values[-1] == 0:
            indices = np.r_[indices, self._sparse_values.indices[-1]]
            values = np.r_[values, 0]
        new_pileup = SparseValues(indices.astype("int"),
                            values.astype("bool"), sanitize=True)
        new_pileup.track_size = self._sparse_values.track_size
        return new_pileup

    def get_node_ids(self):
        node_ids = np.empty_like(self._holes)
        node_ids[:, 0] = np.digitize(self._holes[:, 0], self._node_indexes)
        node_ids[:, 1] = np.digitize(
            self._holes[:, 1], self._node_indexes, True)
        return node_ids

    def _handle_internal(self, internal_holes):
        keep = (internal_holes[:, 1]-internal_holes[:, 0]) > self._max_size
        self._kept_internals = internal_holes[keep]

    def _filter_touched_nodes(self, node_values):
        if not self._touched_nodes:
            return node_values
        touched = [i for i, node in enumerate(node_values[0])
                   if node+self._graph.min_node-1 in self._touched_nodes]
        self._not_touched = np.array(
            [i for i in node_values[0]
             if i+self._graph.min_node-1 not in self._touched_nodes],
            dtype="int")
        return node_values[:, touched]

    def run(self):
        analyzer = SegmentAnalyzer(
            self._holes, self._node_ids, self._node_indexes)
        analyzer.run()
        self._handle_internal(analyzer.internals)
        starts = analyzer.get_ends()
        fulls = self._filter_touched_nodes(analyzer.get_fulls())
        ends = analyzer.get_starts()
        linegraph = DividedLinegraph(
            starts, fulls, ends, self._graph,
            self._last_node, self._touched_nodes)
        mask = linegraph.filter_small(self._max_size)
        self._kept_borders = linegraph.get_masked(mask)
        return self.build_sparse_values(self.build_kept_holes())

    def build_kept_holes(self):
        starts, fulls, ends = self._kept_borders
        fulls = np.r_[fulls, self._not_touched]
        n_starts, n_fulls, n_ends = starts.shape[1], fulls.size, ends.shape[1]
        logging.debug("# %s, %s, %s", n_starts, n_fulls, n_ends)
        n_internals = self._kept_internals.shape[0]
        all_holes = np.empty((n_starts+n_fulls+n_ends+n_internals, 2),
                             dtype="int")
        if all_holes.shape == self._kept_internals.shape:
            return self._kept_internals
        if n_starts:
            all_holes[:n_starts, 0] = self._node_indexes[starts[0]] - starts[1]
            all_holes[:n_starts, 1] = self._node_indexes[starts[0]]
        n_border = n_starts+n_fulls+n_ends
        if fulls.size:
            all_holes[n_starts:n_starts+n_fulls, 0] = self._node_indexes[fulls-1]
            all_holes[n_starts:n_starts+n_fulls, 1] = self._node_indexes[fulls]
        if ends.size:
            all_holes[n_starts+n_fulls:n_border, 0] = self._node_indexes[ends[0]-1]
            all_holes[n_starts+n_fulls:n_border, 1] = self._node_indexes[ends[0]-1] + ends[1]
        all_holes[n_border:] = self._kept_internals
        all_holes.sort(axis=0)
        return all_holes
