import logging
from itertools import chain
import numpy as np
from scipy.stats import poisson
from collections import defaultdict
from .pileup import Pileup
from .pileupcleaner import PileupCleaner
from .pileupcleaner2 import PeaksCleaner, HolesCleaner
from .subgraphcollection import SubgraphCollection
import offsetbasedgraph as obg


class ValuedIndexes(object):
    def __init__(self, indexes, values, start_value, length):
        assert type(indexes) == np.ndarray
        assert type(values) == np.ndarray
        # assert indexes.size == values.size
        self.values = values
        self.indexes = indexes
        self.start_value = start_value
        self.length = length
        self.__tmp_starts = []
        self.__tmp_values = []
        self.__tmp_end = 0

    def __eq__(self, other):
        if np.any(self.values != other.values):
            return False
        if np.any(self.indexes != other.indexes):
            return False

        if isinstance(self.start_value, np.ndarray):
            if np.any(self.start_value != other.start_value):
                return False
        else:
            if self.start_value != other.start_value:
                return False
        return self.length == other.length

    def __str__(self):
        return "VI(%s, %s, %s, %s)" % (
            self.indexes, self.values,
            self.start_value, self.length)

    __repr__ = __str__

    def sum(self):
        lengths = np.diff(self.all_idxs())
        return np.sum(lengths*self.all_values())

    def mean(self):
        return self.sum()/self.length

    def get_subset(self, start, end):
        assert start >= 0
        assert end <= self.length
        if not self.indexes.size:
            return ValuedIndexes(np.array([], dtype="int"),
                                 np.array([]), self.start_value,
                                 end-start)
        length = end-start
        indexes = self.indexes-start
        is_in_subset = np.logical_and(indexes > 0, indexes < length)
        subset_range = np.nonzero(is_in_subset)[0]
        if not subset_range.size:
            return self.__class__(np.array([], dtype="int"),
                                  np.array([]),
                                  self.start_value,
                                  length)

        subset_indexes = indexes[subset_range]
        subset_values = self.values[subset_range]

        if subset_range[0] == 0:
            start_value = self.start_value
        else:
            start_value = self.values[subset_range[0] - 1]
        return self.__class__(subset_indexes,
                              subset_values, start_value,
                              length)

    def scale(self, scale):
        self.values = self.values*scale
        self.start_value = self.start_value*scale

    def tmp_set_interval_value(self, start, end, value):
        self.__tmp_starts.append(start)
        self.__tmp_values.append(value)
        self.__tmp_end = max(self.__tmp_end, end)

    def fix_tmp_values(self):
        assert self.__tmp_starts[0] == 0
        self.start_value = self.__tmp_values[0]
        self.length = self.__tmp_end
        if len(self.__tmp_starts) < 2:
            self.values = np.array([])
            self.indexes = np.array([], dtype="int")
        else:
            self.indexes = np.array(self.__tmp_starts[1:])
            self.values = np.array(self.__tmp_values[1:])
        self.__tmp_values = []
        self.__tmp_starts = []
        self.__tmp_end = 0

    def set_interval_value(self, start, end, value):
        if start == 0:
            self.start_value = value
        else:
            idx = np.nonzero(self.indexes == start)
            self.values[idx] = value

        assert end == self.length or np.any(self.indexes == end)

    def threshold_copy(self, cutoff):
        new_vi = self.__class__(np.copy(self.indexes),
                                self.values >= cutoff,
                                self.start_value >= cutoff,
                                self.length)
        new_vi.sanitize()
        return new_vi

    def threshold(self, cutoff):
        self.values = self.values >= cutoff
        self.start_value = self.start_value >= cutoff

    def sanitize(self):
        diffs = np.diff(np.insert(self.values, 0, self.start_value))
        changes = np.where(diffs != 0)[0]
        self.values = self.values[changes]
        self.indexes = self.indexes[changes]

    def all_values(self):
        return np.insert(self.values, 0, self.start_value)

    def all_idxs(self):
        return np.append(np.insert(self.indexes, 0, 0), self.length)

    def find_valued_areas(self, value):
        all_indexes = self.all_idxs()
        values = self.all_values()
        idxs = np.where(values == value)[0]
        starts = all_indexes[idxs]
        ends = all_indexes[idxs+1]
        return list(chain(*zip(starts, ends)))

    @classmethod
    def maximum(cls, vi_a, vi_b):
        a = vi_a.all_idxs()[:-1]*2
        b = vi_b.all_idxs()[:-1]*2+1
        all_idxs = np.concatenate([a, b])
        all_idxs.sort()
        vi_list = [vi_a, vi_b]
        values_list = []
        for i, vi in enumerate(vi_list):
            idxs = np.nonzero((all_idxs % 2) == i)[0]
            all_values = vi.all_values()
            value_diffs = np.diff(all_values)
            values = np.zeros(all_idxs.shape)
            values[idxs[1:]] = value_diffs
            values[idxs[0]] = all_values[0]
            values_list.append(values.cumsum())
        values = np.maximum(values_list[0], values_list[1])
        idxs = all_idxs // 2
        empty_ends = np.nonzero(np.diff(idxs) == 0)[0]
        max_values = np.maximum(values[empty_ends], values[empty_ends+1])
        values[empty_ends+1] = max_values
        values[empty_ends] = max_values
        obj = cls(idxs[1:], values[1:], values[0], vi_a.length)
        obj.sanitize()
        return obj

    @classmethod
    def combine(cls, vi_a, vi_b):
        a = vi_a.all_idxs()[:-1]*2
        b = vi_b.all_idxs()[:-1]*2+1
        all_idxs = np.concatenate([a, b])
        all_idxs.sort()
        vi_list = [vi_a, vi_b]
        values_list = []
        for i, vi in enumerate(vi_list):
            idxs = np.nonzero((all_idxs % 2) == i)[0]
            all_values = vi.all_values()
            value_diffs = np.diff(all_values)
            values = np.zeros(all_idxs.shape)
            values[idxs[1:]] = value_diffs
            values[0] = all_values[0]
            values_list.append(values.cumsum())

        values = np.array([values_list[0], values_list[1]])
        idxs = all_idxs // 2
        unique_idxs = np.append(np.nonzero(np.diff(idxs))[0], len(idxs)-1)
        idxs = idxs[unique_idxs]
        values = values[:, unique_idxs]
        obj = cls(idxs[1:], np.transpose(values[:, 1:]), values[:, 0], vi_a.length)
        return obj

    def trunctate(self, min_value):
        self.values = np.maximum(self.values, min_value)
        self.start_value = max(min_value, self.start_value)
        self.sanitize()

    @classmethod
    def empty(cls, length):
        return cls(np.array([], dtype="int"),
                   np.array([]),
                   0,
                   length)

    def __iter__(self):
        return zip(
            chain([0], self.indexes),
            chain(self.indexes, [self.length]),
            chain([self.start_value], self.values)
            )

"""
   |     |
    |      |
     |  |
   |       |
"""


class BinaryIndexes(object):
    def __init__(self, starts, ends, length):
        self.starts = starts
        self.ends = ends
        self.length = length

    def add_interval(self, start, end):
        self.starts.append(start)
        self.ends.append(end)
        self.starts.sort


class SparsePileup(Pileup):
    def __init__(self, graph):
        self.graph = graph
        self.data = {rp: ValuedIndexes.empty(graph.node_size(rp))
                     for rp in self.graph.blocks}

    def __eq__(self, other):
        for node_id, vi in other.data.items():
            if self.data[node_id] != vi:
                return False
        return True

    def sum(self):
        return np.sum([values.sum() for node, values in self.data.items()])

    def mean(self):
        graph_size = sum([self.graph.node_size(b) for b in self.graph.blocks])
        mean = self.sum() / graph_size
        return mean

    def scale(self, scale):
        [vi.scale(scale) for vi in self.data.values()]

    def fill_small_wholes(self, max_size):
        # super().fill_small_wholes(max_size)
        cleaner = HolesCleaner(self, max_size)
        areas = cleaner.run()
        for node_id in areas.areas:
            starts = areas.get_starts(node_id)
            ends = areas.get_ends(node_id)
            for start,  end in zip(starts, ends):
                self.data[node_id].set_interval_value(start, end, True)
        self.sanitize()

    def sanitize(self):
        for valued_indexes in self.data.values():
            valued_indexes.sanitize()

    def find_valued_areas(self, value):
        return {node_id: valued_indexes.find_valued_areas(value)
                for node_id, valued_indexes in self.data.items()}

    def set_valued_intervals(self, node_id, valued_indexes):
        self.data[node_id] = valued_indexes

    @classmethod
    def from_areas_collection(cls, graph, areas_list):
        logging.debug(areas_list)
        starts_dict = defaultdict(list)
        ends_dict = defaultdict(list)
        for areas in areas_list:
            for rp in areas.areas:
                starts_dict[rp].extend(areas.get_starts(rp))
                ends_dict[rp].extend(areas.get_ends(rp))
        starts_dict = {rp: np.array(v) for rp, v in starts_dict.items()}
        ends_dict = {rp: np.array(v) for rp, v in ends_dict.items()}
        return cls.from_starts_and_ends(graph, starts_dict,
                                        ends_dict, dtype=int)

    @classmethod
    def from_valued_areas(cls, graph, valued_areas):
        pileup = cls(graph)
        for rp in graph.blocks:
            indexes, values = starts_and_ends_to_sparse_pileup(
                valued_areas.get_starts_array(rp),
                valued_areas.get_ends_array(rp))
            if not indexes.size:
                continue
            start_value = 0
            length = graph.blocks[rp].length()
            if indexes[0] == 0:
                start_value = values[0]
                indexes = indexes[1:]
                values = values[1:]

            if len(indexes) > 0:
                if indexes[-1] == length:
                    indexes = indexes[:-1]
                    values = values[:-1]

            pileup.data[rp] = ValuedIndexes(
                indexes, values, start_value, length)

        return pileup

    @classmethod
    def from_starts_and_ends(cls, graph, starts, ends, dtype=bool):
        pileup = cls(graph)
        for rp in starts:
            indexes, values = starts_and_ends_to_sparse_pileup(
                starts[rp], ends[rp])
            start_value = dtype(False)
            length = graph.blocks[rp].length()
            if indexes[0] == 0:
                start_value = values[0]
                indexes = indexes[1:]
                values = values[1:]

            if len(indexes) > 0:
                if indexes[-1] == length:
                    indexes = indexes[:-1]
                    values = values[:-1]

            pileup.data[rp] = ValuedIndexes(
                indexes, values, start_value, length)

        return pileup

    @classmethod
    def from_intervals(cls, graph, intervals):
        starts, ends = intervals_to_start_and_ends(graph, intervals)
        return cls.from_starts_and_ends(graph, starts, ends)

    def to_subgraphs(self):
        # Returns a list of areas which each is a subgraph
        collection = SubgraphCollection.from_pileup(self.graph, self)
        return collection

    def __str__(self):
        return "\n".join(
            "%s: %s, %s, %d" % (node_id, valued_indexes.indexes, valued_indexes.values, valued_indexes.start_value)
            for node_id, valued_indexes in self.data.items())

    __repr__ = __str__

    def threshold_copy(self, cutoff):
        new_data = {node_id: vi.threshold_copy(cutoff)
                    for node_id, vi in self.data.items()}

        pileup = self.__class__(self.graph)
        pileup.data = new_data
        return pileup

    def threshold(self, cutoff):
        for valued_indexes in self.data.values():
            valued_indexes.threshold(cutoff)
            valued_indexes.sanitize()

    def add_intervals(self, intervals):
        [self.add_interval(interval)
         for interval in intervals]

    def add_interval(self, interval):
        assert all(region_path in self.graph.blocks for
                   region_path in interval.region_paths)
        for i, region_path in enumerate(interval.region_paths):
            start = 0
            end = self.graph.node_size(region_path)
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths)-1:
                end = interval.end_position.offset
            self.data[region_path].add_interval(
                start, end)

    def add_area(self, region_path, start, end, value):
        """
        NAIVE method to add single area to valued_inexes
        Assumes it comes from a complete set of areas
        """
        vi = self.data[region_path]
        vi.tmp_set_interval_value(start, end, value)

    def fix_tmp_values(self):
        for vi in self.data.values():
            vi.fix_tmp_values()

    @classmethod
    def from_bed_graph(cls, graph, filename):
        pileup = super().from_bed_graph(graph, filename)
        assert isinstance(pileup, cls)
        pileup.fix_tmp_values()
        return pileup

    def add_areas(self, areas):
        for area, intervals in areas.items():
            self.data[area].add_interval(
                intervals[0], intervals[1])

    def set_areas_value(self, areas, value):
        for area, intervals in areas.items():
            for i in range(len(intervals)//2):
                self.data[area].set_interval_value(
                    intervals[i], intervals[i+1], value)

    def set_interval_value(self, interval, value):
        assert all(region_path in self.graph.blocks for
                   region_path in interval.region_paths)
        for i, region_path in enumerate(interval.region_paths):
            start = 0
            end = self.graph.node_size(region_path)
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths)-1:
                end = interval.end_position.offset
            self.data[region_path].set_interval_value(
                start, end, value)

    @classmethod
    def from_bed_file(cls, graph, filename):
        f = open(filename, "r")
        starts = defaultdict(list)
        ends = defaultdict(list)
        for line in f:
            if line.startswith("track"):
                continue

            data = line.split()
            block_id = int(data[0])
            starts[block_id].append(int(data[1]))
            ends[block_id].append(int(data[2]))
        starts = {block_id: np.array(start_list) for block_id, start_list in starts.items() if start_list}
        ends = {block_id: np.array(ends_list) for block_id, ends_list in ends.items() if ends_list}
        return cls.from_starts_and_ends(graph, starts, ends)

    def to_bed_file(self, filename):
        f = open(filename, "w")
        areas = self.find_valued_areas(True)
        for node_id, idxs in areas.items():
            for i in range(len(idxs)//2):
                interval = (node_id, idxs[2*i], idxs[2*i+1])
                f.write("%s\t%s\t%s\t.\t.\t.\n" % interval)
        f.close()
        return filename

    def set_to_false(self):
        for valued_indexes in self.data.values():
            valued_indexes.start_value = False
            valued_indexes.indexes = np.array([], dtype="int")
            valued_indexes.values = np.array([], dtype="bool")

    def remove_small_peaks(self, min_size):
        """
        areas = self.find_valued_areas(True)
        intervals = self.areas_to_intervals(areas, include_partial_stubs=False)
        large_intervals = [interval for interval in intervals
                           if interval.length() >= min_size]
        for i in [interval for interval in intervals
                  if interval.length() < min_size]:
            print(i)
        return self.from_intervals(self.graph, large_intervals)
        """
        cleaner = PeaksCleaner(self, min_size)
        areas = cleaner.run()
        pileup = self.from_areas_collection(self.graph, [areas])
        pileup.threshold(0.5)
        return pileup

    def update_max(self, other):
        for key, valued_indexes in self.data.items():
            self.data[key] = ValuedIndexes.maximum(
                valued_indexes, other.data[key])

    def update_max_value(self, min_value):
        for valued_indexes in self.data.values():
            valued_indexes.trunctate(min_value)

    def to_bed_graph(self, filename):
        f = open(filename, "w")
        for node_id, valued_indexes in self.data.items():
            for t in valued_indexes:
                f.write("%s\t%s\t%s\t%s\n" % ((node_id,) + t))
        f.close()
        self.filename = filename
        self.is_written = True
        return filename


class SparseControlSample(SparsePileup):
    def get_p_dict(self):
        p_value_dict = defaultdict(dict)
        count_dict = defaultdict(int)
        for node_id, valued_indexes in self.data.items():
            start = 0
            val = valued_indexes.start_value
            for start, end, val in valued_indexes:
                if val[1] not in p_value_dict[val[0]]:
                    pre_val = poisson.cdf(val[1], val[0])
                    p_val = 1 - pre_val
                    p_value_dict[val[0]][val[1]] = -np.log10(p_val)
                p = p_value_dict[val[0]][val[1]]
                count_dict[p] += end-start

        self.p_value_dict = p_value_dict
        self.count_dict = count_dict

    def get_p_to_q_values(self):
        p_value_counts = self.count_dict
        p_to_q_values = {}
        sorted_p_values = sorted(p_value_counts.keys(), reverse=True)
        rank = 1
        logN = np.log10(sum(p_value_counts.values()))
        pre_q = None
        for p_value in sorted_p_values:
            value_count = p_value_counts[p_value]
            q_value = p_value + (np.log10(rank) - logN)
            if rank == 1:
                q_value = max(0.0, q_value)
            else:
                q_value = max(0.0, min(pre_q, q_value))
            p_to_q_values[p_value] = q_value
            pre_q = q_value
            rank += value_count
        self.p_to_q_values = p_to_q_values

    def get_q_values(self):
        def translation(x):
            return self.p_to_q_values[
                self.p_value_dict[x[0]][x[1]]]
        # f = np.vectorize(translation)
        for valued_indexes in self.data.values():
            if valued_indexes.values.size:
                valued_indexes.values = np.apply_along_axis(
                    translation, 1, valued_indexes.values)

            valued_indexes.start_value = translation(
                valued_indexes.start_value)
            valued_indexes.sanitize()

    def get_scores(self):
        self.get_p_dict()
        self.get_p_to_q_values()
        self.get_q_values()

    @classmethod
    def from_sparse_control_and_sample(cls, control, sample):
        sparse_pileup = cls(control.graph)
        sparse_pileup.data = {
            node_id: ValuedIndexes.combine(
                control.data[node_id], sample.data[node_id])
            for node_id in control.data}
        return sparse_pileup

    @classmethod
    def from_control_and_sample(cls, control, sample):
        sparse_pileup = cls(control.graph)
        for node_id in sparse_pileup.graph.blocks:
            control_counts = control.get_count_array(node_id)
            sample_counts = sample.get_count_array(node_id)
            control_diffs = control_counts[1:]-control_counts[:-1]
            sample_diffs = sample_counts[1:]-sample_counts[:-1]
            changes = np.where(np.logical_or(control_diffs != 0,
                                             sample_diffs))[0] + 1
            vals = np.column_stack([control_counts[changes],
                                    np.floor(sample_counts[changes])])
            init_value = np.array([control_counts[0], sample_counts[0]])
            valued_indexes = ValuedIndexes(
                 changes, vals, init_value, control_counts.shape[0])
            assert not np.any(vals[:, 0] == 0)
            sparse_pileup.set_valued_intervals(node_id, valued_indexes)
        return sparse_pileup


def intervals_to_start_and_ends(graph, intervals):
    # Returns two dict on rp => positions (start/ends)
    starts = defaultdict(list)
    ends = defaultdict(list)

    for interval in intervals:
        for i, region_path in enumerate(interval.region_paths):
            start = 0
            end = graph.node_size(region_path)
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths)-1:
                end = interval.end_position.offset

            if region_path < 0:
                new_start = graph.node_size(region_path) - end
                new_end = graph.node_size(region_path) - start
                start = new_start
                end = new_end

                region_path = -region_path

            starts[region_path].append(start)
            ends[region_path].append(end)

    for rp in starts:
        starts[rp] = np.array(starts[rp])

    for rp in ends:
        ends[rp] = np.array(ends[rp])

    return starts, ends


def filter_pileup_duplicated_position(positions, values):
    equal_to_previous = (positions[:-1] == positions[1:])
    true_array = np.ones(len(positions))
    true_array[np.where(equal_to_previous)] = False
    return positions[np.where(true_array)], values[np.where(true_array)]


def starts_and_ends_to_sparse_pileup(starts, ends):
    coded_starts = starts * 8 + 5
    coded_ends = ends * 8 + 3

    pileup_encoded_positions = np.concatenate((coded_starts, coded_ends))
    pileup_encoded_positions.sort()
    event_codes = (pileup_encoded_positions % 8) - 4  # Containing -1s and 1s
    pileup_values = np.add.accumulate(event_codes)
    pileup_positions = pileup_encoded_positions // 8
    return filter_pileup_duplicated_position(pileup_positions, pileup_values)

