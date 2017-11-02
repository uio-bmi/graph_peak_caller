import numpy as np

def create_control(graph, snarls, reads, d):
    snarl_graph = SnarlGraph(graph, snarls)
    linear_map = LinearSnarlMap(snarl_graph)
    mapped_reads = (linear_map.map_graph_interval(read)
                    for read in reads)
    extended_reads = (((start+end)/2-d, start+end/2+d) for start, end in
                      mapped_reads)
    linear_pileup = LinearPileup.from_starts_and_ends(extended_reads)
    valued_indexes = linear_pileup.to_valued_indexes()
    graph_pileup = SparsePileup(graph)
    graph_pileup.data = valued_indexes
    return graph_pileup


class EventSorter(self):
    def __init__(self, index_lists, values_list, names=None):
        if names is not None:
            [setattr(name.upper(), i) for i, name in enumerate(names)]

        n_types = len(index_lists)
        coded_index_list = [np.array(index_list)*n_types + r for r in
                            range(n_types)]
        coded_indices = np.ravel(coded_index_list)
        sorted_args = np.argsort(coded_indices)
        coded_indices = coded_indices[sorted_args]
        self.indices = coded_indices//n_types
        self.codes = coded_indices % n_types
        self.values = np.array(values_list, dtype="obj").ravel()[sorted_args]

    def __iter__(self):
        return zip(self.coded_indices, self.codes, self.values)


class UnmappedIndices(object):
    def __init__(self):
        self.indices = []
        self.values = []

    def add_indexvalue(self, index, value):
        self.indices.append(index)
        self.values.append(value)


class LinearPileup(object):
    def __init__(self, snarl_graph):
        self._snarl_graph = snarl_graph

    def from_starts_and_ends(self, startends):
        n = len(startends)
        self.starts = np.empty(n, dtype="int")
        self.ends = np.empty(n, dtype="int")
        for i, startend in enumerate(startends):
            self.starts[i] = startend[0]
            self.ends[i] = startend[1]

    def to_valued_indexes(self):
        event_sorter = self.get_event_sorter()
        vi_dict = self.from_event_sorter()
        return vi_dict

    def get_event_sorter(self, linear_map):
        node_start_values = [node_id for node_id in self._snarl_graph.blocks]
        node_end_values = node_start_values[:]
        node_starts_idxs = [linear_map.get_node_start(node_id)
                            for node_id in node_start_values]
        node_end_idxs = [linear_map.get_node_end(node_id)
                         for node_id in node_end_values]
        idxs = [self.indices, node_starts_idxs, node_end_idxs]
        values = [self.values, node_start_values, node_end_values]
        event_sorter = EventSorter(idxs, values, names=["PILEUP_CHANGE",
                                                        "NODE_START",
                                                        "NODE_END"])
        return event_sorter

    def from_event_sorter(self, event_sorter):
        unmapped_indices = {}
        cur_nodes = set([])
        cur_index = 0
        cur_value = 0
        for index, code, value in event_sorter:
            if code == event_sorter.NODE_START:
                cur_nodes.add(value)
                unmapped_indices[value].add_indexvalue(cur_index, cur_value)
            elif code == event_sorter.PILEUP_CHANGE:
                [unmapped_indices[node_id].add_indexvalue(index, value)
                 for node_id in cur_nodes]
                cur_value = value
                cur_index = index
            elif code == event_sorter.NODE_END:
                cur_nodes.remove(value)
        return unmapped_indices

    def to_graph_pileup(self):
        cur_nodes = []
        for idx, event in self.events.items():
            if isinstance(event, tuple):
                if event[0] > 0:
                    cur_nodes.append(event[1])
                else:
                    cur_nodes.remove(event[1])
            else:
                [node.append((idx, value)) for node in cur_nodes]


class LinearSnarlMap(object):
    def __init__(self, snarl_graph):
        self._snarl_graph = snarl_graph
        self._length = self._snarl_graph.length()
        self._linear_node_starts, self._linear_node_ends = snarl_graph.get_distance_dicts()

    def get_node_start(self, node_id):
        return self._linear_node_starts[node_id]

    def get_node_end(self, node_id):
        return self._linear_node_ends[node_id]

    def get_scale_and_offset(self, node_id):
        linear_length = self._linear_node_ends(node_id) - self._linear_node_starts(node_id)
        node_length = self._snarl_graph.node_size(node_id)
        scale = linear_length/node_length
        offset = self._linear_node_starts(node_id)
        return scale, offset

    def to_graph_pileup(self, unmapped_indices_dict):
        vi_dict = {}
        for node_id, unmapped_indices in unmapped_indices_dict.items():
            scale, offset = self.get_scale_and_offset(node_id)
            new_idxs = (np.array(unmapped_indices.indices)-offset)*scale
            new_idxs[0] = min(0, new_idxs[0])
            vi = ValuedIndexes(new_idxs[1:], np.array(unmapped_indices.values)[1:],
                               unmapped_indices.values[0], self._snarl_graph.node_size(node_id))
            vi_dict[node_id] = vi

    def map_graph_interval(self, interval):
        start_pos = self.graph_position_to_linear(interval.start_position)
        end_pos = self.graph_position_to_linear(interval.end_position)
        return start_pos, end_pos

    def graph_position_to_linear(self, position):
        node_start = self._linear_node_starts[position.region_path_id]
        node_end = self._linear_node_ends[position.region_path_id]
        scale = (node_end-node_start)/self.node_sizes(position.node_id)
        return int(node_start + scale*position.offset)

    @classmethod
    def from_snarl_graph(cls, snarl_graph):
        return cls(snarl_graph)
