import offsetbasedgraph as obg
import numpy as np
from collections import deque, defaultdict
import cProfile


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

    def _update_pileup(self, node_id, node_info, starts, prev_ends):
        sub_array = self._pileup.values(node_id)
        node_size = self._graph.node_size(node_id)
        n_in = 0
        for s in node_info._dist_dict.values():
            if s == 0:
                continue
            if s < node_size:
                sub_array[s] -= 1
            # sub_array[0: min(node_size, s)] += 1
            n_in += 1
        sub_array[0] = n_in-prev_ends
        for start in starts:
            sub_array[start] += 1
            if start < node_size-self._fragment_length:
                sub_array[start+self._fragment_length] -= 1

    def run(self):
        start_node = self._graph.get_first_blocks()[0]
        node_info = NodeInfo({})
        queue = deque([(start_node, node_info)])
        unfinished = {}
        cur_id = 0
        while queue:
            node_id, info = queue.popleft()
            node_size = self._graph.node_size(node_id)
            starts = self._starts.get_node_starts(node_id)
            self._update_pileup(node_id, info, starts)
            remains = self._fragment_length - (node_size-starts)
            last_id = cur_id+len(remains)
            remains = dict(zip(range(cur_id, last_id),
                               [r for r in remains if r > 0]))
            remains.update({node: d-node_size for node, d in
                            info._dist_dict.items()
                            if d > node_size})
            cur_id = last_id
            for next_node in self._graph.adj_list[node_id]:
                # print(node_id, next_node)
                if next_node not in unfinished:
                    unfinished[next_node] = TmpNodeInfo(
                        len(self._graph.reverse_adj_list[-next_node]))
                # print("U", unfinished[next_node])
                finished = unfinished[next_node].update(remains)
                if finished is not None:
                    queue.append((next_node, finished))
                    # print("F", finished)
                    del unfinished[next_node]

    def run_linear(self):
        node_ids = self._graph.get_sorted_node_ids()
        prev_ends = 0
        node_infos = defaultdict(NodeInfo)
        cur_id = 0
        empty = NodeInfo()
        for node_id in node_ids:
            info = node_infos.pop(node_id, empty)
            # del node_infos[node_id]
            node_size = self._graph.node_size(node_id)
            starts = self._starts.get_node_starts(node_id)
            self._update_pileup(node_id, info, starts, prev_ends)
            remains = self._fragment_length - (node_size-starts)
            last_id = cur_id + len(remains)
            remains = dict(zip(range(cur_id, last_id),
                               [r for r in remains if r >= 0]))
            remains.update({node: d-node_size for node, d in
                            info._dist_dict.items()
                            if d >= node_size})
            cur_id = last_id
            prev_ends = len(remains)
            for next_node in self._graph.adj_list[node_id]:
                node_infos[next_node].update(remains)
        self._pileup._values = np.cumsum(self._pileup._values)


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
    # print(edges)
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
