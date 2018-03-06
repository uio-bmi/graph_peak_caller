import numpy as np
import logging
from collections import defaultdict

from ..eventsorter import EventSorter, EventSort
from .snarlmaps import LinearSnarlMap


def create_control(linear_map_name, *args, **kwargs):
    logging.info("Reading linear map from file")
    linear_map = LinearSnarlMap.from_json_files(linear_map_name, kwargs["ob_graph"])
    return create_control_from_objs(linear_map, *args, **kwargs)


def create_control_from_objs(linear_map, reads, extension_sizes,
                             fragment_length, ob_graph=None, touched_nodes=None,
                             use_global_min_value=None):
    """
    :param snarl_graph: Hierarchical snarl graph
    """
    linear_size = linear_map._length
    mapped_reads = linear_map.map_interval_collection(reads)
    average_value = mapped_reads.n_intervals*fragment_length / linear_size
    logging.info(
        "Average control value: %.4f (sum of pileup: %d, linear size: %d)" % (
            average_value, mapped_reads.n_basepairs_covered(), linear_size))

    if use_global_min_value is not None:
        average_value = use_global_min_value
        logging.warning("Using global min value %.5f" % use_global_min_value)

    max_pileup = LinearPileup([0], [average_value])
    logging.info("Extending control reads with extension sizes: %s" %
                 extension_sizes)
    for tmp_extension in extension_sizes:
        logging.info("Extension size: %d" % tmp_extension)
        extension = tmp_extension // 2
        extended_reads = mapped_reads.extend(extension)
        linear_pileup = LinearPileup.create_from_starts_and_ends(
                extended_reads.starts, extended_reads.ends)
        assert isinstance(linear_pileup, LinearPileup)
        linear_pileup /= (extension*2/fragment_length)
        logging.info("Linear pileup created. Doing maximum")
        max_pileup.maximum(linear_pileup)

    logging.info("All extensions done. Grating valued indexes from pileup")

    logging.info("Making sparsepilup from valued indexes")
    graph_pileup = max_pileup.to_dense_pileup(linear_map, touched_nodes=touched_nodes)
    logging.info("Control pileup created")

    return graph_pileup


class UnmappedIndices(object):
    def __init__(self, indices=None, values=None):
        self.indices = [] if indices is None else indices
        self.values = [] if values is None else values

    def __str__(self):
        return "(%s, %s)" % (self.indices, self.values)

    def get_index_array(self):
        return np.array(self.indices)

    def get_values_array(self):
        return np.array(self.values)

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
        i = np.array(self.indices)
        v = np.array(self.values)
        pos = i >= 0
        return "Indices: %s, values: %s" % (i[pos],
                                            v[pos])

    __repr__ = __str__

    @classmethod
    def create_from_starts_and_ends(cls, starts, ends):
        logging.info("Creating linear pileup from starts and ends.")
        es = EventSort([starts, ends], [1, -1])
        return LinearPileup(es.indices, es.values)

    def to_sparse_pileup(self, linear_map, touched_nodes=None, min_value=0):
        logging.info("Getting event sorter")
        event_sorter = self.get_event_sorter(linear_map, touched_nodes)
        logging.info("Getting unmapped indices")
        unmapped_indices = self.from_event_sorter(event_sorter)
        logging.info("Mapping linear map to graph pileup")
        #return linear_map.to_numpy_sparse_pileup(unmapped_indices)
        return linear_map.to_sparse_pileup(unmapped_indices, min_value)

    def to_dense_pileup(self, linear_map, touched_nodes=None):
        logging.info("Getting event sorter")
        event_sorter = self.get_event_sorter(linear_map, touched_nodes)
        logging.info("Getting unmapped indices")
        unmapped_indices = self.from_event_sorter(event_sorter)
        logging.info("Mapping linear map to graph pileup")
        return linear_map.to_dense_pileup(unmapped_indices)


    def to_valued_indexes(self, linear_map, touched_nodes=None):
        logging.info("Getting event sorter")
        event_sorter = self.get_event_sorter(linear_map, touched_nodes)
        logging.info("Getting unmapped indices")
        unmapped_indices = self.from_event_sorter(event_sorter)
        logging.info("Mapping linear map to graph pileup")
        vi_dict = linear_map.to_graph_pileup(unmapped_indices)
        return vi_dict

    def get_event_sorter(self, linear_map, touched_nodes=None):
        node_start_values = [node_id for node_id in (linear_map._graph.blocks if touched_nodes is None else touched_nodes)]
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

    @staticmethod
    def from_event_sorter(event_sorter):
        unmapped_indices = defaultdict(UnmappedIndices)
        cur_nodes = set([])
        cur_index = 0
        cur_value = 0
        for index, code, value in event_sorter:
            if code == event_sorter.NODE_START:
                value = int(value)
                cur_nodes.add(value)
                unmapped_indices[value].add_indexvalue(cur_index, cur_value)
            elif code == event_sorter.PILEUP_CHANGE:
                [unmapped_indices[node_id].add_indexvalue(index, value)
                 for node_id in cur_nodes]
                cur_value = value
                cur_index = index
            elif code == event_sorter.NODE_END:
                cur_nodes.remove(int(value))
            else:
                raise Exception("Coding Error")
        return unmapped_indices

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
