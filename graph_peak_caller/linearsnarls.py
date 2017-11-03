import numpy as np
from .snarls import SnarlGraph
from .sparsepileup import SparsePileup, starts_and_ends_to_sparse_pileup
from .snarlmaps import LinearSnarlMap
from .util import sparse_maximum, sanitize_indices_and_values


def create_control(graph, snarls, reads, extension_sizes):
    snarl_graph = SnarlGraph(graph, snarls)
    linear_map = LinearSnarlMap(snarl_graph)

    linear_size = linear_map._length
    mapped_reads = linear_map.map_interval_collection(reads)
    average_value = mapped_reads.n_basepairs_covered() / linear_size

    max_pileup = LinearPileup([], [], snarl_graph)
    for extension in extension_sizes:
        linear_pileup = mapped_reads.extend(extension)
        max_pileup = max_pileup.max(linear_pileup)

    max_pileup.threshold(average_value)


    valued_indexes = max_pileup.to_valued_indexes()
    graph_pileup = SparsePileup(graph)
    graph_pileup.data = valued_indexes
    return graph_pileup


class EventSorter(object):
    def __init__(self, index_lists, values_list, names=None):
        if names is not None:
            [setattr(name.upper(), i) for i, name in enumerate(names)]

        n_types = len(index_lists)
        coded_index_list = [np.array(index_list)*n_types + r 
                            for r in range(n_types)]

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
    def __init__(self, indices, values):
        self.indices = indices
        self.values = values

    @classmethod
    def create_from_starts_and_ends(cls, starts, ends):
        indices, values = starts_and_ends_to_sparse_pileup(starts, ends)
        indices, values = sanitize_indices_and_values(indices, values)
        return LinearPileup(indices, values)

    def to_valued_indexes(self, linear_map):
        event_sorter = self.get_event_sorter()
        unmapped_indices = self.from_event_sorter(event_sorter)
        vi_dict = linear_map.to_graph_pileup(unmapped_indices)
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

    def maximum(self, other):
        indices, values = sparse_maximum(self.indices, self.values,
                                         other.indices, other.values,
                                         max(self.values[-1], other.values[-1]) + 1)
        self.indices = indices
        self.values = values

    def threshold(self, value):
        self.values = np.maximum(self.values, value)


class LinearIntervalCollection(object):

    def __init__(self, starts, ends):
        self.starts = starts
        self.ends = ends

    def extend(self, extension_size):
        extended_starts = (self.starts + self.ends)/2 - extension_size
        extended_ends = (self.starts + self.ends)/2 + extension_size
        linear_pileup = LinearPileup.create_from_starts_and_ends(
                extended_starts, extended_ends)
        return linear_pileup

    def n_basepairs_covered(self):
        return np.sum(self.ends - self.starts)