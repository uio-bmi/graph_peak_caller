import numpy as np
import logging
from scipy.sparse import csr_matrix

from .segmentanalyzer import SegmentSplitter
from ..peakcollection import Peak
from .graphs import PosDividedLineGraph, SubGraph


class SparseMaxPaths:
    def __init__(self, sparse_values, graph, score_pileup, reference_path=None):
        self._node_indexes = graph.node_indexes
        self._sparse_values = sparse_values
        self._graph = graph
        self._score_pileup = score_pileup
        self._segments = self.get_segments()
        self._node_ids = self.get_node_ids()
        self._analyzer = SegmentSplitter(self._segments, self._node_ids,
                                         graph.node_indexes)

        self.reference_mask = None
        if reference_path is not None:
            logging.info("Creating reference mask")
            nodes_in_linear = reference_path.nodes_in_interval()
            logging.info("N nodes before insertion removal: %d" % len(nodes_in_linear))
            # Remove insertions
            i = 0
            for node in nodes_in_linear:
                is_insertion = False
                for prev in self._graph.reverse_adj_list[-node]:
                    if -prev in nodes_in_linear:
                        # Check if -prev has a next on linear
                        for next in self._graph.adj_list[-prev]:
                            if next != node and next in nodes_in_linear:
                                is_insertion = True
                                break
                        break

                if is_insertion:
                    nodes_in_linear.remove(node)
                    logging.info("Removed insertion node %d" % node)

            logging.info("N nodes after insertion removal: %d" % len(nodes_in_linear))

            nodes_in_linear = np.array(list(nodes_in_linear)) - graph.blocks.node_id_offset
            reference_mask = np.zeros_like(self._graph.blocks._array, dtype=bool)
            reference_mask[nodes_in_linear] = True
            self.reference_mask = reference_mask
            logging.info("Done creating reference mask")
        else:
            logging.info("Not creating reference mask")

    def _handle_internal(self, mask):
        ids = self._analyzer._internal_ids[:, 0][mask]
        positions = self._analyzer._internal_segments[mask]
        offsets = positions-self._node_indexes[ids-1, None]
        ids -= self._graph.min_node-1
        logging.info("Found %s internal peaks", len(ids))
        self.internal_paths = [Peak(offset[0], offset[1], [_id],
                                    graph=self._graph)
                               for offset, _id in zip(offsets, ids)]

        for path in self.internal_paths:
            path.info = (False, False)

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

        if self.reference_mask is not None:
            logging.info("Using reference mask when finding max paths to choose path when ambiguous")
            for node_ids, scores in scored_segments:
                is_reference = self.reference_mask[node_ids.astype(int)]
                scores *= 2
                scores[is_reference] += 1
        else:
            logging.info("Not using reference mask when finding max paths")

        self._handle_internal(self._analyzer.internal_mask)
        linegraph = PosDividedLineGraph(scored_segments[2],
                                        scored_segments[3],
                                        scored_segments[1],
                                        self._graph)

        paths, infos, subgraphs = linegraph.max_paths()
        converted = self._convert_paths(paths, infos)
        small_subgraphs = [
            SubGraph(path.region_paths,
                     csr_matrix(([], ([], [])), shape=(1, 1)))
            for path in self.internal_paths]
        return converted+self.internal_paths, subgraphs+small_subgraphs

    def _convert_paths(self, paths, infos):
        reverse_map = np.concatenate(
            [np.flatnonzero(mask) for mask in
             [self._analyzer.end_mask,
              self._analyzer.full_mask,
              self._analyzer.start_mask,
              ]])
        peaks = [self._convert_path(path, reverse_map)
                 for path in paths]
        for peak, info in zip(peaks, infos):
            peak.info = info
        return peaks

    def _convert_path(self, path, reverse_map):
        idxs = reverse_map[path]
        node_ids = self._analyzer._internal_ids[:, 0]
        node_ids = node_ids[idxs]
        start_offset = self._segments[idxs[0], 0 ] - self._node_indexes[node_ids[0]-1]
        end_offset = self._segments[idxs[-1], 1] - self._node_indexes[node_ids[-1]-1]
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
        pileup_cumsum = np.r_[0, np.cumsum(weighted_values)]
        base_scores = pileup_cumsum[pileup_idxs[:, 1]-1]-pileup_cumsum[
            pileup_idxs[:, 0]-1]
        diffs = self._segments-self._score_pileup.indices[pileup_idxs-1]
        values = self._score_pileup.values[pileup_idxs-1]
        val_diffs = diffs*values
        offsets = val_diffs[:, 1] - val_diffs[:, 0]
        self.scores = base_scores + offsets
