import numpy as np
import logging
from .segmentanalyzer import SegmentSplitter
from .peakcollection import Peak
from .sparseholecleaner import PosDividedLineGraph


class SparseMaxPaths:
    def __init__(self, sparse_values, graph, score_pileup):
        self._node_indexes = graph.node_indexes
        self._sparse_values = sparse_values
        self._graph = graph
        self._score_pileup = score_pileup
        self._segments = self.get_segments()
        self._node_ids = self.get_node_ids()
        self._analyzer = SegmentSplitter(self._segments, self._node_ids,
                                         graph.node_indexes)

    def _handle_internal(self, mask):
        ids = self._analyzer._internal_ids[:, 0][mask]
        positions = self._analyzer._internal_segments[mask]
        offsets = positions-self._node_indexes[ids-1, None]
        ids -= self._graph.min_node-1
        logging.info("Found %s internal peaks", len(ids))
        self.internal_paths = [Peak(offset[0], offset[1], [_id],
                                    graph=self._graph)
                               for offset, _id in zip(offsets, ids)]

    def get_segments(self):
        start_idx = 0
        if self._sparse_values.values[0] != 1:
            start_idx += 1
        indices = self._sparse_values.indices[start_idx:]
        if self._sparse_values.values[-1] == 1:
            indices = np.r_[indices, self._node_indexes[-1]]
        return indices.reshape(indices.size//2, 2)

    def get_node_ids(self):
        node_ids = np.empty_like(self._segments)
        node_ids[:, 0] = np.digitize(self._segments[:, 0], self._node_indexes)
        node_ids[:, 1] = np.digitize(self._segments[:, 1],
                                     self._node_indexes, True)
        return node_ids

    def run(self):
        self._analyzer.run()
        self._segments = self._analyzer.splitted_segments
        self.get_segment_scores()
        scored_segments = [np.vstack((self._analyzer._internal_ids[:, 0][mask],
                                      self.scores[mask]))
                           for mask in [self._analyzer.internal_mask,
                                        self._analyzer.start_mask,
                                        self._analyzer.end_mask,
                                        self._analyzer.full_mask]]

        self._handle_internal(self._analyzer.internal_mask)
        linegraph = PosDividedLineGraph(scored_segments[2],
                                        scored_segments[3],
                                        scored_segments[1],
                                        self._graph)

        paths = linegraph.max_paths()
        return self._convert_paths(paths)+self.internal_paths

    def _convert_paths(self, paths):
        reverse_map = np.concatenate(
            [np.flatnonzero(mask) for mask in
             [self._analyzer.end_mask,
              self._analyzer.full_mask,
              self._analyzer.start_mask,
              ]])
        return [self._convert_path(path, reverse_map)
                for path in paths]

    def _convert_path(self, path, reverse_map):
        idxs = reverse_map[path]
        node_ids = self._analyzer._internal_ids[:, 0]
        node_ids = node_ids[idxs]
        start_offset = self._segments[idxs[0], 0 ] - self._node_indexes[node_ids[0]-1]
        end_offset = self._segments[idxs[-1], 1] - self._node_indexes[node_ids[-1]-1]
        # assert np.all(node_ids >= self._graph.blocks.node_id_offset), (self._graph.blocks.node_id_offset,
        # node_ids[node_ids<self._graph.blocks.node_id_offset])
        return Peak(start_offset, end_offset, list(node_ids + self._graph.min_node-1), self._graph)

    def get_segment_scores(self):
        pileup_idxs = np.empty_like(self._segments)
        pileup_idxs[:, 0] = np.digitize(self._segments[:, 0],
                                        self._score_pileup.indices)
        pileup_idxs[:, 1] = np.digitize(self._segments[:, 1],
                                        self._score_pileup.indices, True)
        weighted_values = np.ediff1d(
            self._score_pileup.indices,
            self._score_pileup.track_size-self._score_pileup.indices[-1])*self._score_pileup.values
        # assert np.all(weighted_values >= 0)
        pileup_cumsum = np.r_[0, np.cumsum(weighted_values)]
        base_scores = pileup_cumsum[pileup_idxs[:, 1]-1]-pileup_cumsum[
            pileup_idxs[:, 0]-1]
        assert np.all(base_scores >= 0)
        diffs = self._segments-self._score_pileup.indices[pileup_idxs-1]
        assert np.all(diffs >= 0)
        values = self._score_pileup.values[pileup_idxs-1]
        # assert np.all(values[:, 1] >= values[:, 0])
        val_diffs = diffs*values
        offsets = val_diffs[:, 1] - val_diffs[:, 0]
        self.scores = base_scores + offsets
        assert np.all(self.scores >= 0)
