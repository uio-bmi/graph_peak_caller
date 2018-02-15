from itertools import chain
from collections import defaultdict
import scipy.sparse.csgraph as csgraph
from scipy.sparse import csr_matrix
import numpy as np


class StubsFilter:
    def __init__(self, starts, fulls, ends, graph):
        self._graph = graph
        self._starts = starts
        self._fulls = fulls
        self._ends = ends

        self._starts_mask = np.ones_like(starts, dtype="bool")
        self._fulls_mask = np.ones_like(fulls, dtype="bool")
        self._ends_mask = np.ones_like(ends, dtype="bool")

        self.filter_start_stubs()
        self.filter_end_stubs()

        self.filtered_ends = self._ends[self._ends_mask]
        self.filtered_starts = self._starts[self._starts_mask]
        self.filtered_fulls = self._fulls[self._fulls_mask]
        print(self.filtered_starts)
        print(self.filtered_ends)
        print(self.filtered_fulls)

        self._pos_to_nodes = set(chain(self._fulls, self._ends))
        self._pos_from_nodes = set(chain(self._starts, self._ends))

        self._full_starts = self.find_sub_starts(self.filtered_fulls)
        self._end_starts = self.find_sub_starts(self.filtered_ends)
        self._full_ends = self.find_sub_ends(self.filtered_fulls)
        self._start_ends = self.find_sub_ends(self.filtered_starts)

    def find_sub_starts(self, nodes):
        return np.array([all(-adj in self._pos_from_nodes
                             for adj in self._graph.reverse_adj_list[node])
                         for node in nodes], dtype="bool")

    def find_sub_ends(self, nodes):
        return np.array([all(adj in self._pos_to_nodes
                             for adj in self._graph.adj_list[node])
                         for node in nodes], dtype="bool")

    def _get_start_filter(self, nodes):
        return np.array([bool(self._graph.reverse_adj_list[-node_id])
                         for node_id in nodes], dtype="bool")

    def _get_ends_filter(self, nodes):
        return np.array([bool(self._graph.adj_list[node_id])
                         for node_id in nodes], dtype="bool")

    def filter_start_stubs(self):
        """ Locate nodes that are start_nodes of graph"""
        self._fulls_mask &= self._get_start_filter(self._fulls)
        self._ends_mask &= self._get_start_filter(self._ends)

    def filter_end_stubs(self):
        self._starts_mask &= self._get_ends_filter(self._starts)
        self._fulls_mask &= self._get_ends_filter(self._fulls)


class LineGraph:
    def __init__(self, starts, full, ends, ob_graph):
        self.filtered = StubsFilter(starts[0], full[0], ends[0])



        self.full_nodes = full[0]
        self.start_nodes = starts[0]
        self.end_nodes = ends[0]

        self.full_size = full[1]
        self.start_size = starts[1]
        self.end_size = ends[1]

        self._all_nodes = np.r_[self.start_nodes,
                                self.full_nodes,
                                self.end_nodes]
        self._all_sizes = np.r_[self.start_size,
                                self.full_size,
                                self.end_size]
        self.ob_graph = ob_graph

        self.n_starts = self.start_nodes.size
        self.n_ends = self.start_nodes.size
        self.n_nodes = self._all_nodes.size

        self.mask = np.ones(self.n_nodes)
        self._matrix = self.make_graph()

    def check_stubs(self):
        self.mask[self.n_starts:] &= self.rev_adj_list[-self._all_nodes[self.n_starts:]]
        self.mask[:-self.n_ends] &= self.adj_list[self._all_nodes[:-self.n_ends]]

    def make_graph(self):
        n_starts = self.start_nodes.size
        n_ends = self.end_nodes.size
        to_nodes_dict = {node: n_starts+i for i, node in
                         enumerate(self._all_nodes[n_starts:])}
        self.end_stub = self._all_nodes.size
        print(self.end_stub)
        from_nodes = []
        to_nodes = []
        sizes = []
        possible_to_nodes = set(list(self._all_nodes[n_starts:]))
        print(possible_to_nodes)
        for i, node_id in enumerate(self._all_nodes[:-n_ends]):
            adj_nodes = self.ob_graph.adj_list[node_id]
            if not adj_nodes:
                continue
            print(i, node_id, adj_nodes)
            next_nodes = [to_nodes_dict[node] for node in
                          adj_nodes if node in possible_to_nodes]
            print(next_nodes)
            if len(adj_nodes) != len(next_nodes):
                next_nodes = [self.end_stub]
            from_nodes.extend([i]*len(next_nodes))
            to_nodes.extend(next_nodes)
            sizes.extend([self._all_sizes[i]]*len(next_nodes))

        to_nodes.extend([self.end_stub]*n_ends)
        from_nodes.extend(range(self.end_stub-n_ends, self.end_stub))
        sizes.extend(list(self.end_size))
        print(from_nodes)
        print(to_nodes)
        print(sizes)
        self.n_starts = n_starts
        return csr_matrix((np.array(sizes), (np.array(from_nodes),
                                             np.array(to_nodes))),
                          [self.end_stub+1, self.end_stub+1])

    def filter_small(self, max_size):
        print("#########")
        print(self._matrix)
        print(self._matrix.shape)
        shortest_paths = csgraph.shortest_path(self._matrix)
        to_dist = np.min(shortest_paths[:self.n_starts], axis=0)
        from_dist = shortest_paths[:, self.end_stub]
        print(to_dist+from_dist)
        return to_dist+from_dist <= max_size


class HolesCleaner:
    def __init__(self, graph, sparse_values, max_size):
        self._graph = graph
        self._sparse_values = sparse_values
        self._holes = self.get_holes()
        self._max_size = max_size
        self._node_indexes = graph.node_indexes
        # self._node_size = np.diff(self._node_indexes)

    def get_holes(self):
        start_idx = 0
        if self._sparse_values.values[0] != 0:
            start_idx += 1
        end_idx = self._sparse_values.indices.size
        if self._sparse_values.values[-1] == 0:
            end_idx -= 1

        return self._sparse_values.indices[start_idx:end_idx].reshape(
            (end_idx-start_idx)//2, 2)

    def get_big_holes(self, internal_holes):
        return (internal_holes[:, 1]-internal_holes[0]) > self._max_size

    def run(self):
        border_holes = self.find_border_holes()
        internal_holes = self.holes[~border_holes, :]
        big_internal = self.get_big_holes(internal_holes)
        border_holes = self.holes[border_holes, :]

    def handle_border_holes(self, holes, node_ids):
        # node_offsets = self._node_indexes[node_ids]
        node_id_diffs = node_ids[:, 1]-node_ids[:, 0]
        full_nodes = (chain(range(i[0]+1, i[1])
                            for i in holes[node_id_diffs > 1]))
        full_nodes = [node for node in full_nodes if
                      self._graph.node_size(node) <= self._max_size]

        starts = np.vstack((holes[:, 0], self._node_indexes[node_ids[:, 1]+1]))
        ends = np.vstack((self._node_indexes[node_ids[:, 1]+1], holes[:, 1]))

    def find_border_holes(self):
        holes = self._holes.copy()
        holes[:, 1] += 1
        borders = np.digitize(self._graph.node_indexes, holes.ravel)
        border_holes = np.unique(borders)
        border_holes = border_holes[border_holes%2==1]//2
        return border_holes


if __name__ == "__main__":
    import offsetbasedgraph as obg

    starts = [[5, 10, 19], [10, 100, 3]]
    fulls = [[20, 30, 40, 17, 21], [10, 20, 30, 1, 1]]
    ends = [[50, 53], [60, 1]]
    tmp = {
        5: [20],
        10: [20, 30, 40],
        20: [50],
        30: [50],
        16: [17],
        21: [22],
        53: [54],
    }
    nodes = {i: obg.Block(10) for i in range(5, 100)}
    graph = obg.GraphWithReversals(nodes, tmp)
    StubsFilter(np.array(starts[0]), np.array(fulls[0]),
                np.array(ends[0]), graph)
    # l = LineGraph(np.array(starts), np.array(fulls), np.array(ends), graph)
    # l.filter_small(20)
