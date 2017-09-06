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

        print(type(self.start_value))
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

    def __update_max(self, other):
        sorted(self.tuples+other.tuples, key=lambda x: x[0])
        return None

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


class SparsePileup(Pileup):
    def __init__(self, graph):
        self.graph = graph
        self.data = {}

    def sanitize(self):
        for valued_indexes in self.data.values():
            valued_indexes.sanitize()

    def find_valued_areas(self, value):
        return {node_id: valued_indexes.find_valued_areas(value)
                for node_id, valued_indexes in self.data.items()}

    def set_valued_intervals(self, node_id, valued_indexes):
        self.data[node_id] = valued_indexes

    def __str__(self):
        return "\n".join(
            "%s: %s" % (node_id, sum(array))
            for node_id, array in self.__count_arrays.items())

    __repr__ = __str__

    def set_data(self, count_arrays):
        self.__index_lists = {node_id: [] for node_id in self.graph.blocks}
        self.__value_lists = {node_id: [] for node_id in self.graph.blocks}

    def threshold(self, cutoff):
        for valued_indexes in self.data.values():
            valued_indexes.threshold(cutoff)
            valued_indexes.sanitize()

    def update_max(self, other):
        for node_id, count_array in self.__count_arrays.items():
            self.__count_arrays[node_id] = np.maximum(
                count_array,
                other.__count_arrays[node_id])

    @classmethod
    def from_bed_graph(cls, graph, file_name):
        file = open(file_name)
        pileup = cls(graph)
        # pileup.create_count_arrays()

        graph_uses_int_ids = isinstance(list(graph.blocks.keys())[0], int)

        for line in file:
            if line.startswith("track"):
                continue

            data = line.split()
            block_id = data[0]
            if graph_uses_int_ids:
                block_id = int(block_id)

            start = int(data[1])
            end = int(data[2])
            value = float(data[3])

            pileup.add_area(block_id, start, end, value)

        file.close()
        return pileup

    def create(self):
        self.create_count_arrays()
        for interval in self.intervals:
            self.add_interval(interval)
        return self

    def init_value(self, value):
        self.__count_arrays = {
            node_id: value*np.ones(block.length(), dtype="float")
            for node_id, block in self.graph.blocks.items()}

    def create_count_arrays(self):
        self.__count_arrays = {node_id: np.zeros(block.length(), dtype="float")
                             for node_id, block in self.graph.blocks.items()}

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
            self.valued_intervals[region_path].add_interval(
                start, end)

    def add_area(self, block_id, start, end, value=1):
        count_array_block = self.__count_arrays[block_id]
        self.valued_intervals[block_id].set_interval_value(start, end, value)

    def add_areas(self, areas):
        for area, intervals in areas.items():
            self.valued_intervals[area].add_interval(intervals[0], intervals[1])

    def set_areas_value(self, areas, value):
        for area, intervals in areas.items():
            for i in range(len(intervals)//2):
                self.valued_intervals[area].set_interval_value(
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
            self.valued_intervals[region_path].set_interval_value(
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

    def get_areas_with_value(self):
        areas = {}
        for node_id, count_array in self.__count_arrays.items():
            diffs = count_array[1:]-count_array[:-1]
            changes = np.where(diffs != 0)[0]
            changes = np.pad(changes, 1, "constant")
            changes[0] = -1
            changes += 1
            changes[-1] = self.graph.node_size(node_id)
            vals = count_array[changes[:-1]]
            areas[node_id] = list(zip(changes[:-1], changes[1:], vals))
        return areas

    def to_bed_graph(self, filename):
        f = open(filename, "w")
        areas = self.get_areas_with_value()
        for node_id, tuples in areas.items():
            for t in tuples[:-1]:
                f.write("%s\t%s\t%s\t%s\n" % ((node_id,) + t))
            if tuples[-1][-1] == 0 and not self.graph.adj_list[node_id]:
                continue
            f.write("%s\t%s\t%s\t%s\n" % ((node_id,) + tuples[-1]))
        f.close()
        self.filename = filename
        self.is_written = True
        return filename

    def scale(self, scale):
        for count_array in self.__count_arrays.values():
            count_array *= scale

    def summary(self):
        return sum(array.sum() for array in self.__count_arrays.values())

    def old_find_valued_areas(self, value=False):
        areas = {}
        for node_id, count_array in self.__count_arrays.items():
            cur_area = False
            cur_list = []
            cur_start = None
            for i, val in enumerate(count_array):
                if (val == value) and (not cur_area):
                    cur_start = i
                    cur_area = True
                elif (val != value) and cur_area:
                    cur_list.extend([cur_start, i])
                    cur_area = False
            if cur_area:
                cur_list.extend([cur_start, len(count_array)])
            areas[node_id] = cur_list
        return areas

    def areas_to_intervals(self, areas, include_partial_stubs):
        f = all if include_partial_stubs else any
        intervals = []
        for node_id, startends in areas.items():
            for i in range(len(startends)//2):
                start = startends[2*i]
                end = startends[2*i+1]
                if start == 0:
                    if f(prev_id in areas and areas[prev_id] and areas[prev_id][-1] == self.graph.node_size(prev_id)
                         for prev_id in self.graph.reverse_adj_list[node_id]):
                        continue
                if end == self.graph.node_size(node_id):
                    intervals.extend(
                        self._get_intervals(node_id, [start],
                                            areas, include_partial_stubs))
                else:
                    intervals.append(
                        obg.Interval(start, end, [node_id], graph=self.graph))
        return intervals
            
    def _get_intervals(self, node_id, cur_interval, areas,
                       include_partial_stubs=False):
        intervals = []
        my_interval = cur_interval + [node_id]
        for next_node in self.graph.adj_list[node_id]:
            if (next_node not in areas) or not areas[next_node] or areas[next_node][0] != 0:
                continue
            end = areas[next_node][1]
            if end == self.graph.node_size(next_node):
                intervals.extend(
                    self._get_intervals(next_node, my_interval, areas,
                                        include_partial_stubs))
                continue
            else:
                intervals.append(
                    obg.Interval(
                        my_interval[0], end, my_interval[1:]+[next_node],
                        graph=self.graph))
        f = all if include_partial_stubs else any
        if not f(next_node in areas and areas[next_node] and areas[next_node][0] == 0 for
                 next_node in self.graph.adj_list[node_id]):
            if include_partial_stubs:
                assert self.graph.adj_list[node_id]
            intervals.append(
                obg.Interval(my_interval[0], self.graph.node_size(node_id),
                             my_interval[1:], graph=self.graph))
        return intervals

    def find_valued_intevals(self, value=0):
        interval_dict = {}
        end_interval_dict = {}
        whole_intervals = []
        for node_id, count_array in self.__count_arrays.items():
            cur_area = False
            cur_list = []
            cur_start = None
            for i, val in enumerate(count_array):
                if val == value and not cur_area:
                    cur_start = i
                    cur_area = True
                elif val != value and cur_area:
                    cur_list.append((cur_start, i))
                    cur_area = False
            if cur_area:
                if cur_start == 0:
                    whole_intervals.append(node_id)
                    cur_list = []
                end_interval_dict[node_id] = (cur_start, len(count_array))
            interval_dict[node_id] = cur_list
        return interval_dict, end_interval_dict, whole_intervals

    def fill_small_wholes(self, max_size):
        areas = self.find_valued_areas(False)
        intervals = self.areas_to_intervals(areas, True)
        intervals = [interval for interval in intervals if
                     interval.length() <= max_size]
        for interval in intervals:
            self.set_interval_value(interval, True)

    def remove_small_peaks(self, min_size):
        areas = self.find_valued_areas(True)
        intervals = self.areas_to_intervals(areas, include_partial_stubs=False)
        large_intervals = [interval for interval in intervals
                           if interval.length() >= min_size]
        for node_id, count_array in self.__count_arrays.items():
            count_array *= False
        for interval in large_intervals:
            self.set_interval_value(interval, True)


class SparseControlSample(SparsePileup):
    def get_p_dict(self):
        p_value_dict = {}
        count_dict = defaultdict(int)
        for node_id, valued_indexes in self.data.items():
            start = 0
            val = valued_indexes.start_value
            for start, end, val in valued_indexes:
                sval = str(val)
                if sval not in p_value_dict:
                    p_value_dict[sval] = -np.log10(
                        1 - poisson.cdf(val[1],
                                        val[0]))
                p = p_value_dict[sval]
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
            print(p_value, np.log10(rank), logN, q_value)
            if rank == 1:
                q_value = max(0, q_value)
            else:
                q_value = max(0, min(pre_q, q_value))
            p_to_q_values[p_value] = q_value
            pre_q = q_value
            rank += value_count
        self.p_to_q_values = p_to_q_values

    def get_q_values(self):
        def translation(x):
            return self.p_to_q_values[
                self.p_value_dict[str(x)]]
        # f = np.vectorize(translation)
        for valued_indexes in self.data.values():
            valued_indexes.values = np.apply_along_axis(
                translation, 1, valued_indexes.values)
            valued_indexes.start_value = translation(
                valued_indexes.start_value)

    def get_scores(self):
        print("GETTING P VALUES")
        self.get_p_dict()
        print("GETTING Q VALUES")
        self.get_p_to_q_values()
        print("SETTING Q VALUES")
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
        print("GETTING SPARSE PILEUP")
        sparse_pileup = cls(control.graph)
        for node_id in sparse_pileup.graph.blocks:
            control_counts = control.get_count_array(node_id)
            sample_counts = sample.get_count_array(node_id)
            control_diffs = control_counts[1:]-control_counts[:-1]
            sample_diffs = sample_counts[1:]-sample_counts[:-1]
            changes = np.where(np.logical_or(control_diffs != 0,
                                             sample_diffs))[0] + 1
            vals = np.column_stack([control_counts[changes],
                                    sample_counts[changes]])
            init_value = np.array([control_counts[0], sample_counts[0]])
            valued_indexes = ValuedIndexes(
                 changes, vals, init_value, control_counts.shape[0])
            sparse_pileup.set_valued_intervals(node_id, valued_indexes)
        return sparse_pileup
