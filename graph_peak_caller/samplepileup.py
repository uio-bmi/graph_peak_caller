class StartIndices:
    def __init__(self, indices):
        self._indices = indices


from collections import deque


class NodeInfo:
    def __init__(self, dist_dict):
        self._dist_dict = dist_dict

    def is_ready(self):
        return len(self._dist_dict) == self._edges_in


class TmpNodeInfo(NodeInfo):
    def __init__(self, edges_in):
        self._edges_in = edges_in
        self._dist_dict = {}
        self._n_edges = 0

    def update(self, starts):
        for read_id, val in starts.items():
            self._dist_dict[read_id] = max(self._dist_dict[read_id], val)

        self._n_edges += 1
        if self._n_edges >= self._edges_in:
            return NodeInfo(self._dist_dict)


class PileupCreator:
    def __init__(self, indexed_graph, starts):
        self._graph = indexed_graph
        self._starts = starts
        self._fragment_length = 0

    def _update_pileup(self, node_id, node_info, starts):
        sub_array = self._pileup.get_sub_array(node_id)
        node_size = self._graph.node_size(node_id)
        for s in node_info:
            sub_array[0: min(node_size, s)] += 1
        for start in starts:
            sub_array[start:min(node_size, starts+self._fragment_length)] += 1
                
    def run(self):
        start_node = self._graph.start_node
        node_info = NodeInfo()
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
            remains = dict(zip(range(cur_id, last_id), remains))
            remains.update({node: d-node_size for node, d in
                            info._dist_dict.items()
                            if d > node_size})
            cur_id = last_id
            for next_node in self._graph.adj_list[node_id]:
                if next_node not in unfinished:
                    unfinished[node_id] = TmpNodeInfo(
                        len(self._graph.reverse_adj_list[-next_node]))
                finished = unfinished[node_id].update(remains)
                if finished is not None:
                    queue.append((next_node, finished))
                    del unfinished[next_node]
