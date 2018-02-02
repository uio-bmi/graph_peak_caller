import offsetbasedgraph as obg
import numpy as np
from collections import deque, defaultdict
from itertools import chain
import cProfile
import logging


class StartIndices:
    def __init__(self, indices):
        self._indices = indices


class NodeInfo:
    def __init__(self, dist_dict=None):
        self._dist_dict = dist_dict if dist_dict\
                          is not None else defaultdict(int)

    def update(self, starts):
        for read_id, val in starts.items():
            self._dist_dict[read_id] = max(self._dist_dict[read_id], val)

    def __str__(self):
        return str(self._dist_dict)

    def __repr__(self):
        return str(self._dist_dict)


class TmpNodeInfo(NodeInfo):
    def __init__(self, edges_in):
        self._edges_in = edges_in
        self._dist_dict = defaultdict(int)
        self._n_edges = 0

    def update(self, starts):
        for read_id, val in starts.items():
            self._dist_dict[read_id] = max(self._dist_dict[read_id], val)

        self._n_edges += 1
        if self._n_edges >= self._edges_in:
            return NodeInfo(self._dist_dict)


class PileupCreator:
    def __init__(self, graph, starts, pileup):
        self._graph = graph
        self._starts = starts
        self._fragment_length = 110
        self._pileup = pileup
        self._set_adj_list()

    def _set_adj_list(self):
        self._adj_list = self._graph.adj_list

    def _update_pileup(self, node_info, starts, prev_ends, sub_array):
        node_size = sub_array.size

        n_in = 0
        for s in node_info._dist_dict.values():
            if s == 0:
                continue
            if s < node_size:
                sub_array[s] -= 1
            n_in += 1

        sub_array[0] = n_in-prev_ends
        t_starts = [start for start in starts if start < node_size]
        ends = [start+self._fragment_length for start in starts if
                start < node_size-self._fragment_length]
        for start in t_starts:
            sub_array[start] += 1
        for end in ends:
            sub_array[end] -= 1

    def get_subarray(self, from_idx, to_idx):
        return self._pileup[from_idx:to_idx+1]

    def get_node_ids(self):
        return self._graph.get_sorted_node_ids()

    def run_linear(self):
        node_ids = self.get_node_ids()
        node_infos = defaultdict(NodeInfo)
        cur_id = 0
        empty = NodeInfo()
        cur_array_idx = 0
        i = 0
        for node_id in node_ids:
            if i % 100000 == 0:
                logging.info("%d nodes processed" % i)
            i += 1
            info = node_infos.pop(node_id, empty)
            node_size = self._graph.node_size(node_id)
            starts = self._starts.get_node_starts(node_id)
            n_starts = len(starts)
            n_ends = n_starts + len(info._dist_dict)
            endsidxs = chain(enumerate(
                (start+self._fragment_length for start in starts),
                start=cur_id),
                             info._dist_dict.items())
            self._pileup[cur_array_idx] += n_ends-n_starts
            for start in starts:
                self._pileup[cur_array_idx+start] += 1
            d = {}
            for idx, end in endsidxs:
                if end < node_size:
                    self._pileup[cur_array_idx+end] -= 1
                else:
                    d[idx] = end-node_size
            cur_array_idx += node_size
            cur_id = cur_id + n_starts
            self._pileup[cur_array_idx] -= len(d)
            for next_node in self._adj_list[node_id]:
                node_infos[next_node].update(d)

    def get_pileup(self):
        return np.cumsum(self._pileup)


class ReversePileupCreator(PileupCreator):
    def get_pileup(self):
        return np.cumsum(self._pileup[:-1])[::-1]

    def _set_adj_list(self):
        self._adj_list = self._graph.reverse_adj_list

    def get_node_ids(self):
        return (-node_id for node_id in
                self._graph.get_sorted_node_ids(reverse=True))


class DummyPileup:
    def __init__(self, N, node_size=1000):
        self.N = N
        self._node_size = node_size
        self._values = np.zeros(N*node_size, dtype="int")

    def values(self, k):
        return self._values[(k-1)*node_size:k*node_size]


class DummyStarts:
    def __init__(self, nodes, node_size=1000, dist=100, M=1):
        self._idxs = {node: np.arange(0, node_size, dist)
                      if ((node//100) % M) == 0 else np.array([])
                      for node in nodes}

    def get_node_starts(self, node_id):
        return self._idxs[node_id]


def sorted_wierd_graph(N, node_size):
    nodes = {node_id: obg.Block(node_size)
             for node_id in range(1, 2*N+1)}

    edges = {i: [i-(i % 2)+2, i-(i % 2)+3]
             for i in range(2, 2*N-1)}
    edges[1] = [2, 3]
    return obg.GraphWithReversals(nodes, edges)


def wierd_graph(N, node_size):
    nodes = {node_id: obg.Block(node_size)
             for node_id in range(1, 2*N+1)}
    edges = {i+1: [(i % N)+2, (i % N)+N+2]
             for i in range(2*N-1)}
    edges[N] = []
    edges[2*N] = []
    edges[N+1] = []
    return obg.GraphWithReversals(nodes, edges)

if __name__ == "__main__":
    node_size = 32
    n_nodes = 1000000
    graph = sorted_wierd_graph(n_nodes, node_size)

    # nodes = {node_id: obg.Block(node_size)
    #          for node_id in range(1, n_nodes+1)}
    # edges = {1: [2, 3], 2: [4], 3: [4]}
    # graph = obg.GraphWithReversals(nodes, edges)
    starts = DummyStarts(list(range(1, 2*n_nodes+1)), node_size, dist=10, M=20)
    pileup = DummyPileup(2*n_nodes, node_size)
    creator = PileupCreator(graph, starts, pileup)
    cProfile.run("creator.run_linear()")
