from collections import defaultdict
import logging
import numpy as np

from .sparsediffs import SparseDiffs


def list_array(n):
    array = np.empty(n, dtype=object)
    array[:] = [[] for _ in array]
    return array


class Remains:

    def __init__(self, own_values, node_id):
        own_values = np.sort(own_values)
        self.n_own = own_values.size
        if self.n_own:
            self._remains = {node_id: own_values}
        else:
            self._remains = {}

    def add_remains(self, node_id, remains):
        if node_id in self._remains and self._remains[node_id][-1] > remains[-1]:
            return
        self._remains[node_id] = remains

    def count(self):
        return sum(len(remains) for remains in self._remains.values())

    def add(self, other):
        for node_id, remains in other._remains.items():
            self.add_remains(node_id, remains)

    def __repr__(self):
        return "Remains(%s)" % self._remains

    def cover_distance(self, distance):
        removed = []
        removed_nodes = []
        for key in self._remains:
            cutoff = np.searchsorted(self._remains[key], distance, side="right")
            removed.extend(self._remains[key][:cutoff])
            self._remains[key] = self._remains[key][cutoff:]-distance
            if not self._remains[key].size:
                removed_nodes.append(key)
        for key in removed_nodes:
            del self._remains[key]
        return removed


class EdgePileup:
    def __init__(self, adj_list):
        self._node_index = adj_list._node_index
        self._pileup = np.zeros_like(adj_list._to_nodes)
        self._to_nodes = adj_list._to_nodes

    def add_edge(self, from_node, to_node):
        start, end = self._node_index[from_node:from_node+2]
        self._pileup[start:end] += self._to_nodes[start:end] == to_node

    def add_interval(self, node_ids):
        for from_node, to_node in zip(node_ids[:-1], node_ids[1:]):
            self.add_edge(from_node, to_node)


class MainPileup:
    def __init__(self, graph, linear_repr, d):
        self._node_counter = np.zeros((graph._node_lens.size, 2), dtype="int")
        self._forward_pileup = ForwardPileup(graph, d, self._node_counter)
        self._reverse_pileup = ReversePileup(graph, d, self._node_counter)
        self._edge_pileup = EdgePileup(graph._adj_list)
        self._linear_repr = linear_repr

    def build(self, reads):
        for read in reads:
            self.add_interval(read)
        for pileup in (self._forward_pileup, self._reverse_pileup):
            pileup.extend_intervals()
        return self.combine_strands()

    def linearize_positions(self, position_array):
        linear_pos = [self._linear_repr.get_linear_codes(node_id, offsets)
                      for node_id, offsets in enumerate(position_array)]
        return np.concatenate(linear_pos)

    def combine_strands(self):
        forward_starts = self.linearize_positions(self._forward_pileup.starts)
        reverse_starts = self.linearize_positions(self._reverse_pileup.starts)
        forward_ends = self.linearize_positions(self._forward_pileup.ends)
        reverse_ends = self.linearize_positions(self._reverse_pileup.ends)
        n_starts = len(forward_starts)+len(reverse_starts)
        n_ends = len(forward_ends)+len(reverse_ends)
        all_idxs = np.concatenate(
            ([0],
             forward_starts,
             reverse_starts,
             self._linear_repr.node_offsets[1:-1],
             forward_ends,
             reverse_ends))
        all_diffs = np.empty_like(all_idxs)
        all_diffs[:n_starts+1] = 1
        all_diffs[-n_ends:] = -1
        node_diffs = self._node_counter[1:, 0]-self._node_counter[:-1, 1]
        all_diffs[n_starts+1:-n_ends] = node_diffs
        all_diffs[0] = 0
        args = np.argsort(all_idxs)
        return SparseDiffs(all_idxs[args], all_diffs[args])

    def add_interval(self, interval):
        self._edge_pileup.add_interval(interval.node_ids)
        if interval.direction < 0:
            self._reverse_pileup.add_interval(interval)
        else:
            self._forward_pileup.add_interval(interval)


class Pileup:
    def __init__(self, graph, d, node_counter):
        self.starts = list_array(graph._node_lens.size)
        self.ends = list_array(graph._node_lens.size)
        self._open = list_array(graph._node_lens.size)
        self._graph = graph
        self._d = d

    def extend_intervals(self):
        remains_array = np.array([Remains(ends, node_id)
                                  for node_id, ends in
                                  enumerate(self._open)], dtype=object)

        for node_id, remains in self._items(remains_array):
            remains = remains_array[node_id]
            n_incoming = remains.count()
            self._add_node_starts(node_id, n_incoming-remains.n_own)
            closed = remains.cover_distance(self._graph._node_lens[node_id])
            self._add_node_ends(node_id, n_incoming-len(closed))
            self._set_closed(node_id, closed)
            for next_node in self._adj_list[node_id]:
                remains_array[next_node].add(remains)
            remains_array[node_id] = None

    def _add_node_ids(self, node_ids):
        self._node_starts[node_ids[:-1]] += 1
        self._node_ends[node_ids[1:]] += 1

    def _add_node_starts(self, node_id, count):
        logging.debug("Starts[%s] += %s", node_id, count)
        self._node_starts[node_id] += count

    def _add_node_ends(self, node_id, count):
        logging.debug("Ends[%s] += %s", node_id, count)
        self._node_ends[node_id] += count

    def _add_open(self, node_id, offset):
        self._open[node_id].append(offset)

    def _add_start(self, node_id, offset):
        self.starts[node_id].append(offset)

    def _add_end(self, node_id, offset):
        self.ends[node_id].append(offset)


class ForwardPileup(Pileup):
    def __init__(self, graph, d, node_counter):
        super().__init__(graph, d, node_counter)
        self._node_starts = node_counter[:, 0]
        self._node_ends = node_counter[:, 1]
        self._adj_list = graph._adj_list

    def _set_closed(self, node_id, offsets):
        self.ends[node_id] = offsets

    def add_interval(self, interval):
        self._add_node_ids(interval.node_ids)
        self._add_start(interval.node_ids[0], interval.start)
        self._add_open(interval.node_ids[-1], interval.end+self._d)

    def _items(self, remains_array):
        return enumerate(remains_array)


class ReversePileup(Pileup):
    def __init__(self, graph, d, node_counter):
        super().__init__(graph, d, node_counter)
        self._node_starts = node_counter[:, 1]
        self._node_ends = node_counter[:, 0]
        self._adj_list = self._graph._rev_adj_list

    def add_interval(self, interval):
        self._add_node_ids(interval.node_ids)
        self._add_end(interval.node_ids[-1], interval.end)
        self._add_open(interval.node_ids[0],
                       self._graph._node_lens[interval.node_ids[0]]-interval.start+self._d)

    def _set_closed(self, node_id, offsets):
        self.starts[node_id] = self._graph._node_lens[node_id]-offsets

    def _items(self, remains_array):
        return zip(range(len(remains_array)-1, -1, -1), remains_array[::-1])
