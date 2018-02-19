from .densepileup import DensePileup
from .samplepileup import PileupCreator, ReversePileupCreator
import numpy as np
import logging


class DirectPileup:
    def __init__(self, graph, intervals, pileup):
        self._graph = graph
        self._intervals = intervals
        self._pileup = pileup
        logging.info("Initing index dicts")
        self._pos_ends = {node_id: [] for node_id in self._graph.blocks.keys()}
        self._neg_ends = {-node_id: [] for node_id
                          in self._graph.blocks.keys()}

    def _handle_interval(self, interval):
        self._pileup.add_interval(interval)
        end_pos = interval.end_position
        rp = end_pos.region_path_id
        if rp < 0:
            self._neg_ends[rp].append(end_pos.offset)
        else:
            self._pos_ends[rp].append(end_pos.offset)

    def run(self):
        counter = 0
        for interval in self._intervals:
            self._handle_interval(interval)
            counter += 1

    def to_file(self, base_name):
        indices = np.nonzero(self._pileup.data._values)
        values = self._pileup.data._values[indices]
        np.save(base_name + "_diffindices.npy", indices)
        np.save(base_name + "_diffvalues.npy", values)

    @staticmethod
    def from_file(graph, base_name):
        pileup = DensePileup(graph)
        indices = np.load(base_name + "_diffindices.npy")
        values = np.load(base_name + "_diffvalues.npy")
        pileup.data._values[indices[:-1]] = values[:-1]
        if indices[-1] < pileup.data._values.size:
            pileup.data._values[indices[-1]] = values[-1]
        # pileup.data._values[indices] = values
        pileup.data._values = np.cumsum(pileup.data._values)
        logging.info("Done creating direct pileup from file")
        return pileup


class SparseDirectPileup:
    def add_neg_interval(self, interval):
        rps = [abs(rp)-self.min_id for rp in interval.region_paths]
        starts = self._node_indexes[rps]
        ends = self._node_indexes[1:][rps]
        self._pileup[starts[:-1]] += 1
        self._pileup[ends[1:]] -= 1
        self._pileup[ends[-1]-interval.end_position.offset] += 1
        self._pileup[ends[0]-interval.start_position.offset] -= 1

    def add_interval(self, interval):
        rps = [rp-self.min_id for rp in interval.region_paths]
        starts = self._node_indexes[rps]
        ends = self._node_indexes[1:][rps]
        self._pileup[starts[1:]] += 1
        self._pileup[ends[:-1]] -= 1
        self._pileup[starts[0]+interval.start_position.offset] += 1
        self._pileup[starts[-1] + interval.end_position.offset] -= 1

    def __init__(self, graph, intervals, pileup, out=None):
        self.min_id = pileup.data.min_node
        self._node_indexes = pileup.data._node_indexes
        self._intervals = intervals
        if out is None:
            self._pileup = np.zeros(pileup.data._values.size+1, "int")
        else:
            self._pileup = out
        self._pos_ends = {node_id: [] for node_id in graph.blocks.keys()}
        self._neg_ends = {-node_id: [] for node_id
                          in graph.blocks.keys()}

    def _handle_interval(self, interval):
        # self._pileup.add_interval(interval)
        end_pos = interval.end_position
        rp = end_pos.region_path_id
        if rp < 0:
            self.add_neg_interval(interval)
            self._neg_ends[rp].append(end_pos.offset)
        else:
            self.add_interval(interval)
            self._pos_ends[rp].append(end_pos.offset)

    def run(self):
        i = 0
        for interval in self._intervals:
            if i % 5000 == 0:
                logging.info("%d reads processed" % i)
            self._handle_interval(interval)
            i += 1

    def to_file(self, base_name):
        indices = np.nonzero(self._pileup[:-1])
        values = self._pileup[indices]
        np.save(base_name + "_diffindices.npy", indices)
        np.save(base_name + "_diffvalues.npy", values)

    @staticmethod
    def from_file(graph, base_name):
        pileup = DensePileup(graph)
        indices = np.load(base_name + "_diffindices.npy")
        values = np.load(base_name + "_diffvalues.npy")

        pileup.data._values[indices[:-1]] = values[:-1]
        if indices[-1] < pileup.data._values.size:
            pileup.data._values[indices[-1]] = values[-1]
        pileup.data._values = np.cumsum(pileup.data._values[:-1])
        return pileup


class Starts:
    def __init__(self, d):
        self._dict = d

    def get_node_starts(self, node_id):
        return self._dict[node_id]


# def check_touched2(data):
#     diffs = np.cumsum(data._values)[data._node_indexes[1:]-1]
#     r = np.nonzero(diffs)[0]-1
#     print(data._values[320:400])
#     print(data._node_indexes[10:20])
#     print(diffs[10:20])
#     print(r[10:20])
#     return r

def check_touched(pileup, node_ids):
    return {node_id for node_id in node_ids if
            np.count_nonzero(pileup.data.values(node_id))}


def main(intervals, graph, extension_size, savedirectname=None):
    pileup = DensePileup(graph)
    my_pileup = np.zeros(pileup.data._values.size+2, dtype="int")
    direct_pileup = SparseDirectPileup(graph, intervals, pileup,
                                       out=my_pileup[1:])
    direct_pileup.run()
    if savedirectname is not None:
        direct_pileup.to_file(savedirectname)
    pileup_neg = np.zeros(pileup.data._values.size+1, dtype="int")

    creator = ReversePileupCreator(
        graph, Starts(direct_pileup._neg_ends),
        pileup_neg)
    creator._fragment_length = extension_size
    creator.run_linear()

    my_pileup[1:-1] -= creator._pileup[:0:-1]
    extension_pileup = my_pileup[1:]
    creator = PileupCreator(
        graph, Starts(direct_pileup._pos_ends), extension_pileup)
    creator._fragment_length = extension_size
    creator.run_linear()
    # my_pileup += creator._pileup[:-1]
    pileup.data._values = np.cumsum(my_pileup[1:-1])
    pileup.data._touched_nodes = check_touched(pileup, graph.blocks.keys())
    return pileup
