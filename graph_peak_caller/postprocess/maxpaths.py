import numpy as np
import logging
from scipy.sparse import csr_matrix

from .segmentanalyzer import SegmentSplitter
from ..peakcollection import Peak
from .graphs import PosDividedLineGraph, SubGraph
from .indel_scores import find_invalid_insertions
from .reference_based_max_path import max_path_func


class SparseMaxPaths:
    def __init__(self, sparse_values, graph, score_pileup, reference_path=None, variant_maps=None):
        self._node_indexes = graph.node_indexes
        self._sparse_values = sparse_values
        self._graph = graph
        self._score_pileup = score_pileup
        self._segments = self.get_segments()
        self._node_ids = self.get_node_ids()
        self._analyzer = SegmentSplitter(self._segments, self._node_ids,
                                         graph.node_indexes)

        if variant_maps is not None:
            self._variant_maps = variant_maps
        # if reference_path is not None and variant_maps is not None:
        #     self._invalid_insertions = find_invalid_insertions(score_pileup, variant_maps.insertions, graph, reference_path)
        #     self._insertion_map = variant_maps.insertions
        #     self._reference_path = reference_path
        #     
        # else:
        #     self._invalid_insertions = None
        # self.reference_mask = None
        # if reference_path is not None:
        #     logging.info("Creating reference mask")
        #     nodes_in_linear = reference_path.nodes_in_interval()
        #     logging.info("N nodes before insertion removal: %d" % len(nodes_in_linear))
        #     # Remove insertions
        #     for node in nodes_in_linear:
        #         is_insertion = False
        #         for prev in self._graph.reverse_adj_list[-node]:
        #             if -prev in nodes_in_linear:
        #                 # Check if -prev has a next on linear
        #                 for next in self._graph.adj_list[-prev]:
        #                     if next != node and next in nodes_in_linear:
        #                         is_insertion = True
        #                         break
        #                 break
        # 
        #         if is_insertion:
        #             nodes_in_linear.remove(node)
        #             logging.info("Removed insertion node %d" % node)
        # 
        #     logging.info("N nodes after insertion removal: %d" % len(nodes_in_linear))
        #     
        #     nodes_in_linear = np.array(list(nodes_in_linear)) - graph.blocks.node_id_offset
        #     reference_mask = np.zeros_like(self._graph.blocks._array, dtype=bool)
        #     reference_mask[nodes_in_linear] = True
        #     self.reference_mask = reference_mask
        #     logging.info("Done creating reference mask")
        # else:
        #     logging.info("Not creating reference mask")

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

        # if self.reference_mask is not None:
        #     logging.info("Using reference mask when finding max paths to choose path when ambiguous")
        #     for node_ids, scores in scored_segments:
        #         is_reference = self.reference_mask[node_ids.astype(int)]
        #         scores *= 2
        #         scores[is_reference] += 1
        # else:
        #     logging.info("Not using reference mask when finding max paths")

        if self._invalid_insertions is not None:
            for node_ids, scores in scored_segments:
                print("------------")
                print(self._invalid_insertions[node_ids.astype(int)-1])
                print(scores[self._invalid_insertions[node_ids.astype(int)-1]])
                scores[self._invalid_insertions[node_ids.astype(int)-1]] = 0

        self._handle_internal(self._analyzer.internal_mask)
        linegraph = PosDividedLineGraph(scored_segments[2],
                                        scored_segments[3],
                                        scored_segments[1],
                                        self._graph)
        components = linegraph.get_connected_components()
        components = self._convert_connected_components(components)
        get_max_path = max_path_func(self._score_pileup, self._graph, self._variant_maps)
        max_paths = []
        for i, component in enumerate(components.values()):
            if i % 100 == 0:
                print("path: ", i)
            max_paths.append(get_max_path(component))
        start_nodes = np.array([max_path[0] for max_path in max_paths])-self._graph.min_node
        end_nodes = np.array([max_path[-1] for max_path in max_paths])-self._graph.min_node

        start_args = np.digitize(self._graph.node_indexes[start_nodes+1],
                                 self._sparse_values.indices, right=True)-1
        start_offset = np.maximum(0, self._sparse_values.indices[start_args]-self._graph.node_indexes[start_nodes])
        # prev_indexes = self._sparse_values.indices[start_args-1]
        # node_indexes = self._graph.node_indexes[start_nodes]
        # start_offset = np.where(prev_indexes >= node_indexes, prev_indexes-node_indexes, 0)

        end_args = np.digitize(self._graph.node_indexes[end_nodes], self._sparse_values.indices)-1
        next_indexes = self._sparse_values.indices[end_args+1]
        end_offset = np.minimum(next_indexes, self._graph.node_indexes[end_nodes+1])-self._graph.node_indexes[end_nodes]        
        # max_paths = [get_max_paths(component) for component in components.values()]
        s = SubGraph([], csr_matrix(([], ([], [])), shape=(1, 1)))
        peaks = [Peak(start, end, path, graph=self._graph) for path, start, end in
                 zip(max_paths, start_offset, end_offset)]
        return peaks, [s]*len(max_paths)

        for i in range(10):
            print(components[i])
            print("PATH:", get_max_path(components[i]))
            
        exit()


        paths, infos, subgraphs = linegraph.max_paths()
        for path in paths:
            for node_id in path:
                if self._insertion_map[node_id-1]:
                    print(node_id-1)
        converted = self._convert_paths(paths, infos)
        small_subgraphs = [
            SubGraph(path.region_paths,
                     csr_matrix(([], ([], [])), shape=(1, 1)))
            for path in self.internal_paths]
        return converted+self.internal_paths, subgraphs+small_subgraphs

    def _get_reverse_map(self):
        return np.concatenate(
            [np.flatnonzero(mask) for mask in
             [self._analyzer.end_mask,
              self._analyzer.full_mask,
              self._analyzer.start_mask,
              ]])

    def _convert_paths(self, paths, infos):
        reverse_map = self._get_reverse_map()
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

    def _convert_connected_components(self, node_dict):
        reverse_map= self._get_reverse_map()
        return {key: self._convert_node_ids(node_ids, reverse_map)
                for key, node_ids in node_dict.items()}


    def _convert_node_ids(self,raw_node_ids, reverse_map):
        idxs = reverse_map[raw_node_ids]
        node_ids = self._analyzer._internal_ids[:, 0]
        node_ids = node_ids[idxs]
        return node_ids

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
