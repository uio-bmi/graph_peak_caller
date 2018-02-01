from .densepileup import DensePileup
from .samplepileup import PileupCreator, ReversePileupCreator
import numpy as np


class DirectPileup:
    def __init__(self, graph, intervals, pileup):
        self._graph = graph
        self._intervals = intervals
        self._pileup = pileup
        self._pos_ends = {node_id: [] for node_id in self._graph.blocks.keys()}
        self._neg_ends = {-node_id: [] for node_id in self._graph.blocks.keys()}

    def _handle_interval(self, interval):
        self._pileup.add_interval(interval)
        end_pos = interval.end_position
        rp = end_pos.region_path_id
        if rp < 0:
            self._neg_ends[rp].append(end_pos.offset)
        else:
            self._pos_ends[rp].append(end_pos.offset)

    def run(self):
        for interval in self._intervals:
            self._handle_interval(interval)


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


def main(intervals, graph, extension_size):
    pileup = DensePileup(graph)
    direct_pileup = DirectPileup(graph, intervals, pileup)
    direct_pileup.run()
    pileup_neg = np.zeros_like(pileup.data._values)
    creator = ReversePileupCreator(
        graph, Starts(direct_pileup._neg_ends),
        pileup_neg)
    creator._fragment_length = extension_size
    creator.run_linear()
    pileup.data._values += creator._pileup[::-1]
    del pileup_neg
    extension_pileup = np.zeros_like(pileup.data._values)
    creator = PileupCreator(
        graph, Starts(direct_pileup._pos_ends), extension_pileup)
    creator._fragment_length = extension_size
    creator.run_linear()
    pileup.data._values += creator._pileup
    pileup.data._touched_nodes = check_touched(pileup, graph.blocks.keys())
    return pileup
