from itertools import chain
import numpy as np
import offsetbasedgraph as obg
from collections import defaultdict
from .pileup import Pileup
from scipy.stats import poisson


class ValuedIndexes(object):
    def __init__(self, indexes, values, start_value, length):
        self.values = values
        self.indexes = indexes
        self.start_value = start_value
        self.length = length

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

    def set_interval_value(self, start, end, value):
        if start == 0:
            self.start_value = value
        else:
            idx = np.nonzero(self.indexes == start)
            self.values[idx] = value

        assert end == self.length or np.any(self.indexes == end)

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
        self.data = {}

    def fill_small_wholes(self, max_size):
        super().fill_small_wholes(max_size)
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
    def from_intervals(cls, graph, intervals):
        starts, ends = intervals_to_start_and_ends(graph, intervals)

        pileup = cls(graph)
        for rp in starts:
            indexes, values = starts_and_ends_to_sparse_pileup(starts[rp], ends[rp])
            start_value = False
            length = graph.blocks[rp].length()
            if indexes[0] == 0:
                start_value = values[0]
                indexes = indexes[1:]
                values = values[1:]

            if len(indexes) > 0:
                if indexes[-1] == length:
                    indexes = indexes[:-1]
                    values = values[:-1]

            pileup.data[rp] = ValuedIndexes(indexes, values, start_value, length)

        return pileup

    def __str__(self):
        return "\n".join(
            "%s: %s" % (node_id, valued_indexes.indexes)
            for node_id, valued_indexes in self.data.items())

    __repr__ = __str__

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
        areas = self.find_valued_areas(True)
        intervals = self.areas_to_intervals(areas, include_partial_stubs=False)
        large_intervals = [interval for interval in intervals
                           if interval.length() >= min_size]
        return self.from_intervals(self.graph, large_intervals)


class SparseControlSample(SparsePileup):
    def get_p_dict(self):
        p_value_dict = defaultdict(dict)
        count_dict = defaultdict(int)
        for node_id, valued_indexes in self.data.items():
            start = 0
            val = valued_indexes.start_value
            for start, end, val in valued_indexes:
                # sval = str(val)
                if val[1] not in p_value_dict[val[0]]:
                    p_value_dict[val[0]][val[1]] = -np.log10(
                        1 - poisson.cdf(val[1],
                                        val[0]))
                p = p_value_dict[val[0]][val[1]]
                count_dict[p] += end-start
        # for key, v in p_value_dict.items():
        #    print("%s: %s" % (key, v.keys()))

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

    def to_bed_graph(self, filename):
        f = open(filename, "w")
        for node_id, valued_indexes in self.data.items():
            for t in valued_indexes:
                f.write("%s\t%s\t%s\t%s\n" % ((node_id,) + t))
        f.close()
        self.filename = filename
        self.is_written = True
        return filename

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
    #print("Sorted positions")
    #print(pileup_encoded_positions)

    event_codes = (pileup_encoded_positions % 8) - 4  # Containing -1s and 1s
    #print("Decoded positions")
    #print(event_codes)

    pileup_values = np.add.accumulate(event_codes)
    #print("Accumulation")
    #print(pileup_values)
    #print("Pileup positions")
    pileup_positions = pileup_encoded_positions // 8
    #print(pileup_positions)

    return filter_pileup_duplicated_position(pileup_positions, pileup_values)

