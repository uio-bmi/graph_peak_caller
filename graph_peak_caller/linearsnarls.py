import numpy as np
from collections import defaultdict
from .snarls import SnarlGraph
from .sparsepileup import SparsePileup, starts_and_ends_to_sparse_pileup
from .snarlmaps import LinearSnarlMap
from .util import sparse_maximum, sanitize_indices_and_values
from .eventsorter import EventSorter


def create_control(graph, snarl_graph, reads, extension_sizes):
    """
    :param snarl_graph: Hierarchical snarl graph
    """
    linear_map = LinearSnarlMap(snarl_graph)

    linear_size = linear_map._length
    mapped_reads = linear_map.map_interval_collection(reads)
    average_value = mapped_reads.n_basepairs_covered() / linear_size

    max_pileup = LinearPileup([], [], snarl_graph)
    for extension in extension_sizes:
        extended_reads = mapped_reads.extend()
        linear_pileup = LinearPileup.create_from_starts_and_ends(
                extended_reads.starts, extended_reads.ends)
        max_pileup = max_pileup.max(linear_pileup)

    max_pileup.threshold(average_value)
    valued_indexes = max_pileup.to_valued_indexes()
    graph_pileup = SparsePileup(graph)
    graph_pileup.data = valued_indexes
    return graph_pileup


class UnmappedIndices(object):
    def __init__(self, indices=None, values=None):
        self.indices = [] if indices is None else indices
        self.values = [] if values is None else values

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
        event_sorter = self.get_event_sorter(linear_map)
        unmapped_indices = self.from_event_sorter(event_sorter)
        vi_dict = linear_map.to_graph_pileup(unmapped_indices)
        return vi_dict

    def get_event_sorter(self, linear_map):
        node_start_values = [node_id for node_id in linear_map._graph.blocks]
        node_end_values = node_start_values[:]
        node_starts_idxs = [linear_map.get_node_start(node_id)
                            for node_id in node_start_values]
        node_end_idxs = [linear_map.get_node_end(node_id)
                         for node_id in node_end_values]
        idxs = [node_end_idxs, self.indices, node_starts_idxs]
        values = [node_end_values, self.values, node_start_values]
        event_sorter = EventSorter(idxs, values, names=["NODE_END",
                                                        "PILEUP_CHANGE",
                                                        "NODE_START",
                                                        ])
        return event_sorter

    def from_event_sorter(self, event_sorter):
        unmapped_indices = defaultdict(UnmappedIndices)
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
