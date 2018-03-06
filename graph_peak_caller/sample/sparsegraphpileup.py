import numpy as np
import logging
from itertools import chain
from collections import defaultdict

import offsetbasedgraph as obg

from ..sparsediffs import SparseDiffs
from ..densepileup import DensePileup


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


class SparseGraphPileup:
    def __init__(self, graph):
        self.starts = []
        self.ends = []
        self.node_starts = np.zeros(len(graph.blocks)+2)
        self.touched_nodes = np.zeros(self.node_starts.size, dtype="bool")

    def __str__(self):
        return("\n".join([str(self.starts), str(self.ends), str(self.node_starts)]))


class ReadsAdder:
    def __init__(self, graph, pileup, keep_ends=False):
        self.min_id = graph.min_node
        self._node_indexes = graph.node_indexes
        self._pileup = pileup
        self._pos_ends = {node_id: [] for node_id in graph.blocks.keys()}
        self._neg_ends = {-node_id: [] for node_id
                          in graph.blocks.keys()}

    def _handle_interval(self, interval):
        # self._pileup.add_interval(interval)
        end_pos = interval.end_position
        rp = end_pos.region_path_id
        if rp < 0:
            self.add_open_neg_interval(interval)
            self._neg_ends[rp].append(end_pos.offset)
        else:
            self.add_open_pos_interval(interval)
            self._pos_ends[rp].append(end_pos.offset)

    def get_pos_ends(self):
        return self._pos_ends

    def get_neg_ends(self):
        return self._neg_ends

    def _add_start_pos(self, pos):
        rp = (abs(pos.region_path_id)-self.min_id)
        if pos.region_path_id > 0:
            self._pileup.starts.append(self._node_indexes[rp]+pos.offset)
        else:
            self._pileup.ends.append(self._node_indexes[rp+1]-pos.offset)

    def add_open_pos_interval(self, interval):
        rps = [rp-self.min_id for rp in interval.region_paths]
        self._pileup.touched_nodes[rps] = True
        self._pileup.node_starts[rps[1:]] += 1
        self._pileup.node_starts[1:][rps[:-1]] -= 1
        self._add_start_pos(interval.start_position)

    def add_open_neg_interval(self, interval):
        rps = [abs(rp)-self.min_id for rp in interval.region_paths]
        self._pileup.touched_nodes[rps] = True
        self._pileup.node_starts[1:][rps[1:]] -= 1
        self._pileup.node_starts[rps[:-1]] += 1
        self._add_start_pos(interval.start_position)

    def add_reads(self, intervals):
        i = 0
        for interval in intervals:
            if i % 5000 == 0:
                logging.info("%d reads processed" % i)
            self._handle_interval(interval)
            i += 1


class ReadsAdderWDirect(ReadsAdder):
    def __init__(self, graph, pileup):
        self.pos_read_ends = []
        self.neg_read_ends = []
        super().__init__(graph, pileup)

    def _handle_interval(self, interval):
        # self._pileup.add_interval(interval)
        super()._handle_interval(interval)
        end_pos = interval.end_position
        rp = end_pos.region_path_id
        if rp < 0:
            self.neg_read_ends.append(
                self._node_indexes[abs(rp)+1-self.min_id]-end_pos.offset)
        else:
            self.pos_read_ends.append(
                self._node_indexes[abs(rp)-self.min_id]+end_pos.offset)


class SparseExtender:
    def __init__(self, graph, pileup, fragment_length):
        self._graph = graph
        self._pileup = pileup
        self._fragment_length = fragment_length
        self._set_adj_list()
        self._graph_size = graph.node_indexes[-1]

    def _set_adj_list(self):
        self._adj_list = self._graph.adj_list

    def get_node_ids(self):
        return range(self._graph.min_node, self._graph._id+1)
    # return self._graph.get_sorted_node_ids()

    def _add_node_end(self, node_id, value):
        self._pileup.node_starts[node_id+1-self._graph.min_node] -= value

    def _add_node_start(self, node_id, value):
        self._pileup.node_starts[node_id-self._graph.min_node] += value

    def _add_end(self, index):
        self._pileup.ends.append(index)

    def run_linear(self, starts_dict):
        node_ids = self.get_node_ids()
        node_infos = defaultdict(NodeInfo)
        cur_id = 0
        empty = NodeInfo()
        cur_array_idx = 0
        n_nodes = len(node_ids)
        for i, node_id in enumerate(node_ids):
            if i % 1000000 == 0:
                logging.info("Handling node %s of %s", i, n_nodes)
            info = node_infos.pop(node_id, empty)
            if info._dist_dict:
                self._pileup.touched_nodes[
                    abs(node_id)-self._graph.min_node] = True
            node_size = self._graph.node_size(node_id)
            starts = starts_dict[node_id]
            n_starts = len(starts)
            # n_ends = n_starts + len(info._dist_dict)
            endsidxs = chain(enumerate(
                (start+self._fragment_length for start in starts),
                start=cur_id),
                             info._dist_dict.items())
            self._add_node_start(node_id, len(info._dist_dict))
            d = {}
            for idx, end in endsidxs:
                if end <= node_size:
                    self._add_end(cur_array_idx+end)
                else:
                    d[idx] = end-node_size
            cur_array_idx += node_size
            cur_id = cur_id + n_starts
            self._add_node_end(node_id, len(d))
            for next_node in self._adj_list[node_id]:
                node_infos[next_node].update(d)

    def get_pileup(self):
        return np.cumsum(self._pileup)


class ReverseSparseExtender(SparseExtender):
    def get_pileup(self):
        return np.cumsum(self._pileup[:-1])[::-1]

    def _set_adj_list(self):
        self._adj_list = self._graph.reverse_adj_list

    def get_node_ids(self):
        return range(-self._graph._id, -self._graph.min_node+1)
    # return (-node_id for node_id in
    #             self._graph.get_sorted_node_ids(reverse=True))

    def _add_node_end(self, node_id, value):
        self._pileup.node_starts[-node_id-self._graph.min_node] += value

    def _add_node_start(self, node_id, value):
        self._pileup.node_starts[-node_id+1-self._graph.min_node] -= value

    def _add_end(self, index):
        self._pileup.starts.append(self._graph_size-index)


def _get_node_indexes(graph):
    if isinstance(graph.blocks, obg.BlockArray):
        # Quicker way to make node_indexes array
        logging.info("(using cumsum on np block array)")
        node_indexes = np.cumsum(graph.blocks._array, dtype=np.uint32)
        logging.info("Node indexes created...")
        graph.min_node = (graph.blocks.node_id_offset+1)
        return node_indexes
    sorted_nodes = sorted(graph.blocks.keys())
    min_node = sorted_nodes[0]
    graph.min_node = min_node
    max_node = sorted_nodes[-1]
    span = max_node-min_node+1
    node_indexes = np.zeros(span+1, dtype=np.uint32)
    offset = 0
    for i, node in enumerate(sorted_nodes):
        index = node - min_node
        node_indexes[index] = offset
        offset += graph.node_size(node)
        node_indexes[-1] = offset
    return node_indexes


class SamplePileupGenerator:
    def __init__(self, graph, extension):
        graph.node_indexes = _get_node_indexes(graph)
        self._pileup = SparseGraphPileup(graph)
        self._graph = graph
        self._reads_adder = ReadsAdderWDirect(graph, self._pileup)
        self._pos_extender = SparseExtender(graph, self._pileup, extension)
        self._neg_extender = ReverseSparseExtender(
            graph, self._pileup, extension)

    def __to_dense(self):
        diffs = SparseDiffs.from_pileup(self._pileup,
                                        self._graph.node_indexes)
        sparse_values = diffs.get_sparse_values()
        densepileup = DensePileup(self._graph)
        for start, end, value in zip(
                sparse_values.indices,
                chain(sparse_values.indices[1:], [densepileup.data._values.size]),
                sparse_values.values):
            densepileup.data._values[start:end] = value
        densepileup.data._touched_nodes = {
            node_id for node_id in self._graph.blocks if
            np.count_nonzero(densepileup.data.values(node_id))}
        return densepileup

    def save_direct(self, base_name):
        new_pileup = SparseGraphPileup(self._graph)
        new_pileup.ends = self._pileup.ends + self._reads_adder.pos_read_ends
        new_pileup.starts = self._pileup.starts + self._reads_adder.neg_read_ends
        new_pileup.node_starts = self._pileup.node_starts
        sparsediff = SparseDiffs.from_pileup(
            new_pileup, self._graph.node_indexes)
        sparsediff.clean()
        sparse_values = sparsediff.get_sparse_values()
        sparse_values.track_size = self._graph.node_indexes[-1]
        sparse_values.to_sparse_files(base_name)
        # np.save(base_name + "_diffindices.npy", sparsediff._indices)
        # np.save(base_name + "_diffvalues.npy", sparsediff._diffs)

    def run(self, reads, save_name=None):
        self._reads_adder.add_reads(reads)
        if save_name is not None:
            self.save_direct(save_name)
        self._pos_extender.run_linear(self._reads_adder.get_pos_ends())
        self._neg_extender.run_linear(self._reads_adder.get_neg_ends())
        sdiffs = SparseDiffs.from_pileup(self._pileup,
                                         self._graph.node_indexes)
        sdiffs.touched_nodes = set(
            np.flatnonzero(
                self._pileup.touched_nodes[:-2]) + self._graph.min_node)
        return sdiffs
    # return self.__to_dense()
