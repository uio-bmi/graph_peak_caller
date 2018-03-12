import numpy as np
import logging
from .linearintervals import LinearIntervalCollection
from ..sparsediffs import SparseDiffs


class LinearMap:
    def __init__(self, node_starts, node_ends, graph):
        self._node_starts = np.asanyarray(node_starts)
        self._node_ends = np.asanyarray(node_ends)
        assert np.all(self._node_starts < self._node_ends),\
            np.where(self._node_starts >= self._node_ends)
        self._length = self._node_ends[-1]+graph.node_size(graph.min_node)
        self._graph = graph

    def __eq__(self, other):
        if not np.all(self._node_starts == other._node_starts):
            return False
        return np.all(self._node_ends == other._node_ends)

    def __repr__(self):
        return "%s:%s" % (self._node_starts, self._node_ends)

    def get_node_start(self, node_id):
        return self._node_starts[node_id-self._graph.min_node]

    def get_node_end(self, node_id):
        return self._node_ends[node_id-self._graph.min_node]

    def get_scale_and_offset(self, node_id):
        linear_length = self.get_node_end(node_id) \
                        - self.get_node_start(node_id)
        node_length = self._graph.node_size(node_id)
        scale = linear_length/node_length
        offset = self.get_node_start(node_id)
        return scale, offset

    def map_interval_collection(self, interval_collection):
        starts = []
        ends = []
        for interval in interval_collection:
            start, end = self.map_graph_interval(interval)
            starts.append(start)
            ends.append(end)
        return LinearIntervalCollection(starts, ends)

    def map_graph_interval(self, interval):
        start_pos = self.graph_position_to_linear(interval.start_position)
        end_pos = self.graph_position_to_linear(interval.end_position)
        return start_pos, end_pos

    def graph_position_to_linear(self, position):
        node_id = abs(position.region_path_id)-self._graph.min_node
        node_start = self._node_starts[node_id]
        node_end = self._node_ends[node_id]
        node_size = self._graph.node_size(position.region_path_id)
        scale = (node_end-node_start) / node_size
        if position.region_path_id > 0:
            return node_start + scale*position.offset
        else:
            return node_end - scale*position.offset

    def to_sparse_pileup(self, unmapped_indices_dict, min_value=0):
        all_indices = []
        all_values = []
        i = 0
        node_idxs = self._graph.node_indexes
        min_idx = self._graph.min_node
        for i in range(node_idxs.size-1):
            if i % 100000 == 0:
                logging.info("Processing node %d" % i)
            node_id = i+min_idx
            if node_id not in unmapped_indices_dict:
                all_indices.append(node_idxs[i])
                all_values.append(min_value)
                continue
            unmapped_indices = unmapped_indices_dict[node_id]
            scale, offset = self.get_scale_and_offset(node_id)
            new_idxs = [(idx-offset)//scale+node_idxs[i]
                        for idx in unmapped_indices.indices]
            new_idxs[0] = max(node_idxs[i], new_idxs[0])
            all_indices.extend(new_idxs)
            all_values.extend(unmapped_indices.values)
        return SparseDiffs(
            np.array(all_indices, dtype="int"),
            np.diff(np.r_[0, all_values]))

    @classmethod
    def from_graph(cls, graph):
        starts = cls.find_starts(graph)
        ends = cls.find_ends(graph)
        return cls(starts, ends, graph)

    @classmethod
    def from_file(cls, filename, graph):
        obj = np.load(filename)
        return cls(obj["starts"], obj["ends"], graph)

    def to_file(self, filename):
        np.savez(filename, starts=self._node_starts,
                 ends=self._node_ends)

    @staticmethod
    def find_starts(graph):
        node_ids = list(graph.get_sorted_node_ids())
        max_dists = np.zeros(len(node_ids))
        counter = 0
        for i, node_id in enumerate(node_ids):
            assert graph.node_size(node_id) > 0
            cur_dist = max_dists[i] + graph.node_size(node_id)
            if (i) and max_dists[i] == 0:
                # print(i, node_id, graph.reverse_adj_list[-node_id])
                # raise
                pass
            for next_node in graph.adj_list[node_id]:
                if next_node < node_id:
                    counter += 1
                j = next_node-graph.min_node
                max_dists[j] = max(cur_dist, max_dists[j])
        print("### Counts", counter)
        print("###", max_dists[-1]+graph.node_size(node_ids[-1]))
        print(np.max(max_dists))
        return max_dists

    @staticmethod
    def find_ends(graph):
        adj_list = graph.reverse_adj_list
        node_ids = list(graph.get_sorted_node_ids(reverse=True))
        max_dists = np.zeros(len(node_ids))
        N = len(node_ids)
        for i, node_id in enumerate(node_ids):
            cur_dist = max_dists[N-i-1] + graph.node_size(node_id)
            assert (not i) or max_dists[N-i-1] > 0, i
            for next_node in adj_list[-node_id]:
                j = abs(next_node)-graph.min_node
                max_dists[j] = max(cur_dist, max_dists[j])
        assert np.argmax(max_dists) == 0
        linear_length = max_dists[0] + graph.node_size(graph.min_node)
        print("###", linear_length)
        return linear_length - max_dists

if __name__ == "__main__":
    import offsetbasedgraph as obg
    graph = obg.GraphWithReversals.from_numpy_file("graph.nobg")
    LinearMap.from_graph(graph)
