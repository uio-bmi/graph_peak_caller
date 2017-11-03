import numpy as np
from .sparsepileup import ValuedIndexes


class LinearSnarlMap(object):
    def __init__(self, snarl_graph, graph):
        self._snarl_graph = snarl_graph
        self._graph = graph
        self._length = self._snarl_graph.length()
        self._linear_node_starts, self._linear_node_ends = snarl_graph.get_distance_dicts()

    def get_node_start(self, node_id):
        return self._linear_node_starts[node_id]

    def get_node_end(self, node_id):
        return self._linear_node_ends[node_id]

    def get_scale_and_offset(self, node_id):
        linear_length = self.get_node_end(node_id) - self.get_node_start(node_id)
        node_length = self._graph.node_size(node_id)
        scale = linear_length/node_length
        offset = self.get_node_start(node_id)
        return scale, offset

    def to_graph_pileup(self, unmapped_indices_dict):
        vi_dict = {}
        for node_id, unmapped_indices in unmapped_indices_dict.items():
            scale, offset = self.get_scale_and_offset(node_id)
            new_idxs = (np.array(unmapped_indices.indices)-offset) * scale
            new_idxs = new_idxs.astype("int")
            new_idxs[0] = min(0, new_idxs[0])
            vi = ValuedIndexes(
                new_idxs[1:], np.array(unmapped_indices.values)[1:],
                unmapped_indices.values[0], self._graph.node_size(node_id))
            vi_dict[node_id] = vi
        return vi_dict

    def map_graph_interval(self, interval):
        start_pos = self.graph_position_to_linear(interval.start_position)
        end_pos = self.graph_position_to_linear(interval.end_position)
        return start_pos, end_pos

    def graph_position_to_linear(self, position):
        node_id = abs(position.region_path_id)
        node_start = self._linear_node_starts[node_id]
        node_end = self._linear_node_ends[node_id]
        node_size = self._graph.node_size(position.region_path_id)
        scale = (node_end-node_start) / node_size
        if position.region_path_id > 0:
            return node_start + scale*position.offset
        else:
            return node_end - scale*position.offset

    @classmethod
    def from_snarl_graph(cls, snarl_graph):
        return cls(snarl_graph)

    def map_interval_collection(self, interval_collection):
        starts = []
        ends = []
        for interval in interval_collection:
            start, end = self.map_graph_interval(interval)
            starts.append(start)
            ends.append(end)
        return LinearIntervalCollection(starts, ends)
