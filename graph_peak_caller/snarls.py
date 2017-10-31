import offsetbasedgraph as obg
from pyvg import Snarls
from collections import deque, defaultdict


class SnarlGraph(obg.GraphWithReversals):

    def __init__(self, blocks, edges, start_node, end_node):
        self._length = None
        self._blocks = blocks
        self._edges = edges
        self._start_node = start_node
        self._end_node = end_node
        super(SnarlGraph, self).__init__(blocks, edges)
        self._get_linear_start_and_end_pos()

    def get_next_nodes(self, node_id):
        if node_id not in self._edges:
            return []
        return self._edges[node_id]

    def get_previous_nodes(self, node_id):
        return self.reverse_adj_list[node_id]

    def node_size(self, node_id):
        return self._blocks[abs(node_id)].length()

    @classmethod
    def create_from_start_and_end(cls, start, end, parent_graph):
        pass

    def length(self):
        if self._length is not None:
            return self._length
        self._length = self._get_longest_path_length()
        return self._length

    def _get_longest_path_length(self):
        memo = self._create_path_length_dict()
        return memo[self._end_node]

    def _create_path_length_dict(self, forward=True):
        start_node = self._start_node if forward else -self._end_node
        end_node = self._end_node if forward else -self._start_node
        next_node_func = self.get_next_nodes if forward else self.get_previous_nodes
        stack = deque([(start_node, 0)])
        memo = defaultdict(int)
        while stack:
            node_id, dist = stack.pop()
            for next_node in next_node_func(node_id):
                if memo[next_node] > dist:
                    continue
                memo[next_node] = dist
                if next_node == end_node:
                    continue
                new_dist = dist + self.node_size(next_node)
                stack.append((next_node, new_dist))

        return memo

    def _get_linear_start_and_end_pos(self):
        self._forward_length_dict = self._create_path_length_dict()
        self._back_length_dict = self._create_path_length_dict(False)

    def get_distance_dicts(self):
        starts_dict, ends_dict = {}, {}
        for node_id, block in self._blocks.items():
            # TODO: Also get reverse distances
            if isinstance(block, obg.Block):
                starts_dict[node_id] = self._forward_length_dict[node_id]
                ends_dict[node_id] = self._back_length_dict[-node_id]
            else:
                sub_starts_dict, sub_ends_dict = block.get_distance_dicts()
                starts_dict.update(
                    {n_id: self._forward_length_dict[node_id]+dist
                     for n_id, dist in sub_starts_dict.items()})
                ends_dict.update({n_id: self._back_length_dict[-node_id]+dist
                                  for n_id, dist in sub_ends_dict.items()})

        return starts_dict, ends_dict


class SimpleSnarl():
    def __init__(self, start, end, id, parent=None):
        self.start = start
        self.end = end
        self.id = id
        self.parent = parent
        

class SnarlGraphBuilder():

    def __init__(self, graph, snarls):
        self.snarls = snarls

    @classmethod
    def from_vg_snarls(cls, vg_snarls_file_name):

        snarls = Snarls.from_vg_snarls_file(vg_snarls_file_name)

        for snarl in snarls:
            if hasattr(snarl, "parent"):
                print(snarl)


