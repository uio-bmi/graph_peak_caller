import numpy as np
from .sparse_pileup import SparsePileup
from offsetbasedgraph import Position


class ShiftParameters(object):
    def __init__(self, mfold, winsize=10):
        self.mfold = mfold
        self.winsize = winsize


class PeakModel(object):
    def __init__(self, graph, intervals, shift_parameters):
        self.graph = graph
        self.intervals = intervals
        self.parameters = shift_parameters
        self.pileups = {}

    def _create_pileups(self):
        for direction in [+1, -1]:
            self.pileups[direction] = SparsePileup.from_intervals(
                (i for i in self.intervals if i.direction == direction))

    def _find_peaks(self):
        for direction in [+1, -1]:
            self.peaks[direction] = self.pileups[direction].threshold(
                self.parameters.mfold)

    def _find_centers(self):
        intervals = {direction: peaks.find_valued_intevals(True)
                     for direction, peaks in self.peaks.items()}

        self.centers = {
            direction: [self.pileups[direction].find_max(interval)
                        for interval in intervals]
            for direction, intervals in intervals.items()}

    def _rec_find_pairs_on_edge(self):
        for node_id, centers_array in self.centers_arrays[node]:
            last_center = centers_array[-1]
            if last_center % 3 != 0:
                continue
            last_center //= 3
            if self.graph.node_size(node_id) - last_center > self.parameters.win_size:
                pass

    def _find_pairs_on_edge(self):
        pairs = []
        for node, edge_list in self.graph.adj_list.items():
            pos_center = self.centers_arrays[node][-1]
            if pos_center % 3 != 0:
                continue
            pos_center //= 3
            distance_to_end = self.graph.node_size(node) - pos_center
            if distance_to_end > self.parameters.win_size:
                continue
            for next_node in edge_list:
                neg_center = self.centers_arrays[node]
                if start % 3 != 1:
                    continue
                neg_center //= 3
                if neg_center + distance_to_end <= self.parameters.win_size:
                    pairs.append(Position(node, pos_center),
                                 Position(next_node, neg_center))
        return pairs

    def _find_paired_peaks_on_node(self, centers_array):
        centers_array.sort()
        diffs = np.diff(centers_array)
        pair_idxs = np.where(diffs % 3 == 1 and diffs < (self.params.win_size * 3))[0]
        return centers_array[pair_idxs], centers_array[pair_idxs+1]

    def _find_paired_peaks(self):
        centers_dicts = {}
        for direction in [+1, -1]:
            centers_dicts[direction] = defaultdict(list)
            for center in self.centers[direction]:
                centers_dicts[direction][center.region_path].append(
                    center.offset)

        centers_arrays = {node_id: np.array(self.centers_dicts[+1][node_id] +
                                            self.centers_dicts[-1][node_id])
                          for node_id in set(centers_dicts[+1].keys() + centers_dicts[-1].keys())}

        for node_id, center_array in centers_arrays.items():
            centers_arrays *= 3
            centers_arrays[len(self.centers_dicts[-1]):] += 1

        pairs = []
        for node_id, array in centers_arrays.items():
            pos_centers, neg_centers = self._find_paired_peaks(array)
            pairs.extend([
                (Position(node_id, pos_center),
                 Position(node_id, neg_center))
                for pos_center, neg_center in
                zip(pos_centers, neg_centers)])

        pairs.extend(self._find_pairs_on_edge())

    def _get_windows(self):
        self.
