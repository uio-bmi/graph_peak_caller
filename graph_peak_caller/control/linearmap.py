import numpy as np


class LinearMap:
    def __init__(self, node_starts, node_ends, graph):
        self._node_starts = node_starts
        self._node_ends = node_ends
        self._graph = graph

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
        # node_sizes = graph.blocks._array
        for i, node_id in enumerate(node_ids):
            cur_dist = max_dists[i] + graph.node_size(node_id)
            for next_node in graph.adj_list[node_id]:
                j = next_node-graph.min_node
                max_dists[j] = max(cur_dist, max_dists[j])
        return max_dists

    @staticmethod
    def find_ends(graph):
        adj_list = graph.reverse_adj_list
        node_ids = list(graph.get_sorted_node_ids(reverse=True))
        max_dists = np.zeros(len(node_ids))
        N = len(node_ids)
        for i, node_id in enumerate(node_ids):
            cur_dist = max_dists[N-i-1] + graph.node_size(node_id)
            for next_node in adj_list[-node_id]:
                j = abs(next_node)-graph.min_node
                max_dists[j] = max(cur_dist, max_dists[j])
        linear_length = max_dists[0]+graph.node_size(graph.min_node)
        return linear_length - max_dists

if __name__ == "__main__":
    import offsetbasedgraph as obg
    graph = obg.GraphWithReversals.from_numpy_file("graph.nobg")
    LinearMap.from_graph(graph)
