import logging
from itertools import chain
import scipy.sparse.csgraph as csgraph
import numpy as np
from scipy.sparse import csr_matrix
from .subgraphanalyzer import SubGraphAnalyzer


class DummyTouched:
    def __contains__(self, item):
        return True


class StubsFilter:
    def __init__(self, starts, fulls, ends, graph,
                 last_node=None, touched_nodes=None):
        self._last_node = last_node
        self._graph = graph
        self._starts = starts
        self._fulls = fulls
        self._ends = ends

        self._touched_nodes = touched_nodes
        if touched_nodes is None:
            self._touched_nodes = DummyTouched()

        self._starts_mask = np.ones_like(starts, dtype="bool")
        self._fulls_mask = np.ones_like(fulls, dtype="bool")
        self._ends_mask = np.ones_like(ends, dtype="bool")

        self.filter_start_stubs()
        self.filter_end_stubs()
        self.filtered_ends = self._ends[self._ends_mask]
        self.filtered_starts = self._starts[self._starts_mask]
        self.filtered_fulls = self._fulls[self._fulls_mask]
        self._set_pos_nodes()

        self._full_starts = np.flatnonzero(
            self.find_sub_starts(self.filtered_fulls))
        self._end_starts = np.flatnonzero(
            self.find_sub_starts(self.filtered_ends))
        self._full_ends = np.flatnonzero(
            self.find_sub_ends(self.filtered_fulls))
        self._start_ends = np.flatnonzero(
            self.find_sub_ends(self.filtered_starts))

    def find_sub_starts(self, nodes):
        r = np.array([not all(-adj in self._pos_from_nodes or adj
                              not in self._touched_nodes
                              for adj in self._graph.reverse_adj_list[-node])
                      for node in nodes], dtype="bool")
        return r

    def find_sub_ends(self, nodes):
        a = np.array([not all(adj in self._pos_to_nodes or adj > self._last_node or adj not in self._touched_nodes
                              for adj in self._graph.adj_list[node])
                      for node in nodes], dtype="bool")
        return a

    def _get_start_filter(self, nodes):
        return np.array([len(self._graph.reverse_adj_list[-node_id])
                         for node_id in nodes], dtype="bool")

    def _get_ends_filter(self, nodes):
        return np.array([len(self._graph.adj_list[node_id])
                         for node_id in nodes], dtype="bool")

    def filter_start_stubs(self):
        """ Locate nodes that are start_nodes of graph"""
        self._fulls_mask &= self._get_start_filter(self._fulls)
        self._ends_mask &= self._get_start_filter(self._ends)

    def filter_end_stubs(self):
        self._starts_mask &= self._get_ends_filter(self._starts)
        self._fulls_mask &= self._get_ends_filter(self._fulls)

    def _set_pos_nodes(self):
        self._pos_to_nodes = set(chain(self._fulls, self._ends))
        self._pos_from_nodes = set(chain(self._starts, self._fulls))


class PosStubFilter(StubsFilter):

    def _set_pos_nodes(self):
        self._pos_to_nodes = set(chain(self._fulls, self._ends))
        self._pos_from_nodes = set(chain(self._starts, self._fulls))

    def find_sub_starts(self, nodes):
        return np.array([not any(-adj in self._pos_from_nodes
                                 for adj in self._graph.reverse_adj_list[-int(node)])
                         for node in nodes], dtype="bool")

    def find_sub_ends(self, nodes):
        return np.array([not any(adj in self._pos_to_nodes
                                 for adj in self._graph.adj_list[int(node)])
                         for node in nodes], dtype="bool")

    def filter_start_stubs(self):
        """ Locate nodes that are start_nodes of graph"""
        pass

    def filter_end_stubs(self):
        pass


class SubGraph:
    def __init__(self, node_ids, graph):
        self._node_ids = node_ids
        self._graph = graph

    def __str__(self):
        return "Subgraph (%s)" % self._node_ids

    def __repr__(self):
        return self.__str__()


class LineGraph:
    stub_class = StubsFilter
    def __init__(self, starts, full, ends, ob_graph, last_node=None, touched_nodes=None):
        if last_node is not None:
            last_node += ob_graph.min_node-1
        starts[0] += ob_graph.min_node-1
        full[0] += ob_graph.min_node-1
        ends[0] += ob_graph.min_node-1
        self.filtered = self.stub_class(starts[0], full[0], ends[0],
                                        ob_graph, last_node, touched_nodes)
        self.full_size = full[1][self.filtered._fulls_mask]
        self.start_size = starts[1][self.filtered._starts_mask]
        self.end_size = ends[1][self.filtered._ends_mask]

        self.kept_starts = np.vstack((starts[0][~self.filtered._starts_mask],
                                      starts[1][~self.filtered._starts_mask]))
        self.kept_ends = np.vstack((ends[0][~self.filtered._ends_mask],
                                    ends[1][~self.filtered._ends_mask]))
        self.kept_fulls = full[0][~self.filtered._fulls_mask]
        self.kept_fulls = self.kept_fulls.astype("int")

        self.full_nodes = self.filtered.filtered_fulls
        self.start_nodes = self.filtered.filtered_starts
        self.end_nodes = self.filtered.filtered_ends
        self._all_nodes = np.r_[self.start_nodes,
                                self.full_nodes,
                                self.end_nodes]
        self._all_nodes = self._all_nodes.astype("int")
        self._all_sizes = np.r_[self.start_size,
                                self.full_size,
                                self.end_size]
        self._all_sizes = self._all_sizes.astype("int")
        self.ob_graph = ob_graph

        self.n_starts = self.start_nodes.size
        self.n_ends = self.end_nodes.size
        self.n_nodes = self._all_nodes.size
        self._matrix = self.make_graph()

    def make_graph(self):
        n_starts = self.start_nodes.size
        n_ends = self.end_nodes.size
        to_nodes_dict = {node: n_starts+i for i, node in
                         enumerate(self._all_nodes[n_starts:])}
        self.end_stub = self._all_nodes.size
        from_nodes = []
        to_nodes = []
        sizes = []
        possible_to_nodes = set(list(self._all_nodes[n_starts:]))
        for i, node_id in enumerate(self._all_nodes[:self._all_nodes.size-n_ends]):
            adj_nodes = self.ob_graph.adj_list[node_id]
            next_nodes = [to_nodes_dict[node] for node in
                          adj_nodes if node in possible_to_nodes]
            from_nodes.extend([i]*len(next_nodes))
            to_nodes.extend(next_nodes)
            sizes.extend([self._all_sizes[i]]*len(next_nodes))

        end_nodes = np.r_[self.filtered._start_ends,
                          n_starts + self.filtered._full_ends,
                          np.arange(self.n_nodes-n_ends, self.n_nodes)]
        to_nodes.extend([self.end_stub]*end_nodes.size)
        from_nodes.extend(list(end_nodes))
        sizes.extend(self._all_sizes[end_nodes])
        # self.n_starts = n_starts
        self._graph_node_sizes = sizes
        return csr_matrix((np.array(sizes, dtype="int32"),
                           (np.array(from_nodes, dtype="int32"),
                            np.array(to_nodes, dtype="int32"))),
                          [self.end_stub+1, self.end_stub+1])

    def get_masked(self, mask):
        n_starts, n_ends = self.n_starts, self.n_ends
        self._all_nodes -= (self.ob_graph.min_node-1)
        T = self._all_nodes.size
        starts = np.vstack((self._all_nodes[:n_starts][mask[:n_starts]],
                            self._all_sizes[:n_starts][mask[:n_starts]]))
        fulls = self._all_nodes[n_starts:T-n_ends][mask[n_starts:T-n_ends]]
        ends = np.vstack((self._all_nodes[T-n_ends:][mask[T-n_ends:]],
                          self._all_sizes[T-n_ends:][mask[T-n_ends:]]))
        d = self.ob_graph.min_node-1
        self.kept_starts[0] -= d
        self.kept_ends[0] -= d
        self.kept_fulls -= d
        starts = np.hstack((starts, self.kept_starts))
        ends = np.hstack((ends, self.kept_ends))
        fulls = np.r_[fulls, self.kept_fulls]

        return starts, fulls, ends

    def filter_small(self, max_size):
        start_nodes = np.r_[np.arange(self.n_starts),
                            self.n_starts+self.filtered._full_starts,
                            self.n_nodes-self.n_ends+self.filtered._end_starts]
        if not start_nodes.size:
            return np.array([], dtype="bool")
        connected_components = csgraph.connected_components(self._matrix)
        shortest_paths = csgraph.shortest_path(self._matrix)
        to_dist = np.min(shortest_paths[start_nodes], axis=0)
        from_dist = shortest_paths[:, self.end_stub]
        return (to_dist+from_dist)[:-1] > max_size


class DividedLinegraph(LineGraph):

    def _get_subgraph(self, idxs):
        self._lookup[idxs] = np.arange(idxs.size)
        sub_matrix = self._matrix[idxs]
        new_indices = self._lookup[sub_matrix.indices]
        sub_matrix = csr_matrix((sub_matrix.data,
                                 new_indices,
                                 sub_matrix.indptr),
                                shape=(idxs.size, idxs.size))
        return sub_matrix

    def filter_small(self, max_size):
        start_nodes = np.r_[np.arange(self.n_starts),
                            self.n_starts+self.filtered._full_starts,
                            self.n_nodes-self.n_ends+self.filtered._end_starts]
        if not start_nodes.size:
            return np.ones_like(self._all_nodes, dtype="bool")
        start_nodes_mask = np.zeros(self.end_stub, dtype="bool")
        start_nodes_mask[start_nodes] = True
        if not start_nodes.size:
            return np.array([], dtype="bool")
        n_components, connected_components = csgraph.connected_components(
            self._matrix[:self.end_stub, :self.end_stub])
        logging.info("Found %s components", n_components)
        complete_mask = np.zeros(self.end_stub, dtype="bool")
        self._lookup = np.empty(self.end_stub+1, dtype="int")
        for comp in range(n_components):
            if comp % 1000 == 0:
                logging.info("Component %s of %s", comp, n_components)

            idxs = np.r_[np.flatnonzero(
                connected_components == comp), self.end_stub]

            if idxs.size-1 > 36:
                complete_mask[idxs[:-1]] = True
                continue
            my_mask = start_nodes_mask[idxs[:-1]]
            start_nodes = np.flatnonzero(my_mask)
            if not start_nodes.size:
                complete_mask[idxs[:-1]] = True
                continue
            subgraph = self._get_subgraph(idxs)
            shortest_paths = csgraph.shortest_path(subgraph)
            to_dist = np.min(shortest_paths[start_nodes], axis=0)
            from_dist = shortest_paths[:, -1]
            complete_mask[idxs[:-1]] = (to_dist+from_dist)[:-1] > max_size
        return complete_mask


class PosDividedLineGraph(DividedLinegraph):
    stub_class = PosStubFilter

    def _backtrace(self, dists, predecessors, start_nodes):
        start = start_nodes[np.argmin(dists[start_nodes, -1])]
        max_row = predecessors[start]
        cur = max_row[-1]
        path = []
        while cur > 0:
            path.append(cur)
            cur = max_row[cur]
        if (not path) or path[-1] != start:
            path += [start]
        return path

    def max_paths(self):
        start_nodes = np.r_[np.arange(self.n_starts),
                            self.n_starts+self.filtered._full_starts,
                            self.n_nodes-self.n_ends+self.filtered._end_starts]
        if not start_nodes.size:
            return [], [], []
        start_nodes_mask = np.zeros(self.end_stub, dtype="bool")
        start_nodes_mask[start_nodes] = True
        paths = []
        infos = []
        self._matrix.data += 1
        subgraphs = []
        n_components, connected_components = csgraph.connected_components(
            self._matrix[:self.end_stub, :self.end_stub])
        logging.info("Found %s components", n_components)
        # node_id_refs = []
        # subgraphs = []
        for comp in range(n_components):
            if comp % 100 == 0:
                logging.info("Component %s of %s", comp, n_components)
            idxs = np.r_[np.flatnonzero(
                connected_components == comp), self.end_stub]
            subgraph = -self._matrix[idxs][:, idxs]
            subgraph.data += 1
            subgraphs.append(SubGraph(self._all_nodes[idxs[:-1]], subgraph))
            if idxs.size == 2:
                paths.append([idxs[0]])
                infos.append((0, 0))
                continue
            distances, predecessors = csgraph.shortest_path(
                subgraph, return_predecessors=True, method="BF")
            my_mask = start_nodes_mask[idxs[:-1]]
            start_nodes = np.flatnonzero(my_mask)
            local_idxs = self._backtrace(distances, predecessors, start_nodes)
            global_idxs = idxs[local_idxs]
            paths.append(global_idxs[::-1])
            infos.append(
                SubGraphAnalyzer(subgraph, local_idxs[::-1], distances).get_info())
        return paths, infos, subgraphs

