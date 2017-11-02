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
        self._get_linear_mapped_node_intervals()

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
        print(self._forward_length_dict)
        print(self._back_length_dict)

    def _get_linear_mapped_node_intervals(self):
        self.linear_node_intervals = {}
        for node_id in self._blocks:
            start_length = self._forward_length_dict[node_id]
            path_length = start_length + self._back_length_dict[-node_id] + self.node_size(node_id)
            scale_factor = self.length()/path_length
            print(node_id, path_length, start_length, scale_factor)
            self.linear_node_intervals[node_id] = (
                start_length * scale_factor,
                (start_length+self.node_size(node_id))*scale_factor)

    def get_distance_dicts(self):
        starts_dict, ends_dict = {}, {}
        for node_id, block in self._blocks.items():
            # TODO: Also get reverse distances
            start, end = self.linear_node_intervals[node_id]
            if isinstance(block, obg.Block):
                print("B: ", node_id, start, end)
                starts_dict[node_id] = start
                ends_dict[node_id] = end
            else:
                sub_starts_dict, sub_ends_dict = block.get_distance_dicts()
                print("#", node_id, start)
                print(sub_ends_dict)
                starts_dict.update(
                    {n_id: start+dist
                     for n_id, dist in sub_starts_dict.items()})
                ends_dict.update({n_id: start + dist
                                  for n_id, dist in sub_ends_dict.items()})

        return starts_dict, ends_dict

    def find_node_ids(self, offset):
        sub_nodes = [node_id for node_id, interval
                     in self.linear_node_intervals.items()
                     if interval[0] < offset and interval[1] > offset]
        positions = []
        for node_id in sub_nodes:
            if isinstance(self._blocks[node_id],
                          obg.Block):
                positions.append(node_id)
            else:
                positions.extend(self._blocks[node_id].find_node_ids())

        return positions


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


