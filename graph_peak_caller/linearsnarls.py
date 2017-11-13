import numpy as np
from collections import defaultdict
from .sparsepileup import SparsePileup
from .util import sanitize_indices_and_values
from .eventsorter import EventSorter, EventSort
from .snarlmaps import LinearSnarlMap
import logging


def create_control(linear_map_name, reads, extension_sizes, fragment_length):
    """
    :param snarl_graph: Hierarchical snarl graph
    """
    linear_map = LinearSnarlMap.from_file(linear_map_name)
    linear_size = linear_map._length
    mapped_reads = linear_map.map_interval_collection(reads)
    average_value = mapped_reads.n_intervals*fragment_length / linear_size
    logging.info("Average control value: %.4f (sum of pileup: %d, linear size: %d)" % (average_value, mapped_reads.n_basepairs_covered(), linear_size))
    max_pileup = LinearPileup([0], [average_value])

    for tmp_extension in extension_sizes:
        extension = tmp_extension//2
        print("#########", extension)
        extended_reads = mapped_reads.extend(extension)
        linear_pileup = LinearPileup.create_from_starts_and_ends(
                extended_reads.starts, extended_reads.ends)
        linear_pileup /= (extension*2/fragment_length)
        print(linear_pileup)
        max_pileup.maximum(linear_pileup)
    # max_pileup.threshold(average_value)
    valued_indexes = max_pileup.to_valued_indexes(linear_map)
    graph_pileup = SparsePileup(linear_map._graph)
    graph_pileup.data = valued_indexes

    return graph_pileup


class UnmappedIndices(object):
    def __init__(self, indices=None, values=None):
        self.indices = [] if indices is None else indices
        self.values = [] if values is None else values

    def __str__(self):
        return "(%s, %s)" % (self.indices, self.values)

    def add_indexvalue(self, index, value):
        self.indices.append(index)
        self.values.append(value)


class LinearPileup(object):
    def __init__(self, indices, values):
        self.indices = indices
        self.values = values

    def __eq__(self, other):
        if not np.allclose(self.indices, other.indices):
            return False
        return np.allclose(self.values, other.values)

    def __itruediv__(self, scalar):
        self.values /= scalar
        return self

    def __str__(self):
        pos = self.indices >= 0
        return "Indices: %s, values: %s" % (self.indices[pos],
                                            self.values[pos])

    __repr__ = __str__

    @classmethod
    def create_from_starts_and_ends(cls, starts, ends):
        es = EventSort([starts, ends], [1, -1])
        return LinearPileup(es.indices, es.values)

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
        for start_idx, end_idx in zip(node_starts_idxs, node_end_idxs):
            assert start_idx < end_idx
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
                try:
                    cur_nodes.remove(value)
                except:
                    raise
            else:
                raise Exception("Coding Error")
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
                [node.append((idx, value))
                 for node in cur_nodes]

    def continuous_sparse_maximum(self, other):
        indices1 = self.indices
        indices2 = other.indices
        values1 = self.values
        values2 = other.values
        all_indices = np.concatenate([indices1, indices2])
        codes = np.concatenate([np.zeros_like(indices1),
                                np.ones_like(indices2)])
        sorted_args = np.argsort(all_indices)
        sorted_indices = all_indices[sorted_args]
        sorted_codes = codes[sorted_args]
        values_list = []
        for code, values in enumerate((values1, values2)):
            my_args = np.where(sorted_codes == code)[0]
            diffs = np.diff(values)
            my_values = np.zeros(sorted_indices.shape)
            my_values[my_args[1:]] = diffs
            my_values[my_args[0]] = values[0]
            values_list.append(my_values.cumsum())
        values = np.maximum(values_list[0], values_list[1])
        self.indices = sorted_indices
        self.values = values
        self.sanitize_indices()
        self.sanitize_values()
        # 
        # 
        # empty_ends = np.nonzero(np.diff(sorted_indices) == 0)[0]
        # max_values = np.maximum(values[empty_ends], values[empty_ends+1])
        # values[empty_ends+1] = max_values
        # values[empty_ends] = max_values
        # indices, values = sanitize_indices_and_values(sorted_indices, values)
        # self.indices = indices
        # self.values = values

    def sanitize_indices(self, choose_last=True):
        assert choose_last
        idx_diffs = np.diff(self.indices)
        changes = np.nonzero(idx_diffs)[0]
        new_args = np.concatenate([changes, [self.values.size - 1]])
        self.indices = self.indices[new_args]
        self.values = self.values[new_args]

    def sanitize_values(self):
        value_diffs = np.diff(self.values)
        changes = np.nonzero(value_diffs)[0]
        new_args = np.concatenate([[0], changes+1])
        self.indices = self.indices[new_args]
        self.values = self.values[new_args]

    def maximum(self, other):
        return self.continuous_sparse_maximum(other)

    def threshold(self, value):
        self.values = np.maximum(self.values, value)
        self.sanitize_values()
