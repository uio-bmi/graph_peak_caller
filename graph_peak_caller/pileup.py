from itertools import chain
import numpy as np
import offsetbasedgraph as obg


class Pileup(object):
    def __init__(self, graph):
        self.graph = graph
        self.filename = None
        self.is_written = False
        if self.graph is not None:
            self.create_count_arrays()

    def __str__(self):
        return "\n".join(
            "%s: %s" % (node_id, sum(array))
            for node_id, array in self.__count_arrays.items())

    __repr__ = __str__

    def set_count_arrays(self, count_arrays):
        self.__count_arrays = count_arrays

    def __eq__(self, other):
        if False and self.graph != other.graph:
            return False
        for key, value in self.__count_arrays.items():
            if key not in other.__count_arrays:
                return False
            if not all(value == other.__count_arrays[key]):
                return False

        return len(self.__count_arrays) == len(other.__count_arrays)

    def is_numerically_equal(self, other):
        if False and self.graph != other.graph:
            return False
        for key, value in self.__count_arrays.items():
            if key not in other.__count_arrays:
                return False
            if not all(np.isclose(value, other.__count_arrays[key])):
                return False

        return len(self.__count_arrays) == len(other.__count_arrays)

    def threshold(self, cutoff):
        for node_id, count_array in self.__count_arrays.items():
            self.__count_arrays[node_id] = count_array > cutoff
            assert self.__count_arrays[node_id].dtype == np.dtype("bool")

    def clear(self):
        self.__count_arrays = None

    def get_count_arrays(self):
        return self.__count_arrays

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
        [self.add_interval(interval) for interval in intervals]

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
            self.__count_arrays[region_path][start:end] += 1

    def add_area(self, block_id, start, end, value=1):
        count_array_block = self.__count_arrays[block_id]
        count_array_block[start:end] += value

    def add_areas(self, areas):
        for area, intervals in areas.items():
            self.__count_arrays[area][intervals[0]:intervals[1]] += 1
            # for i in range(len(intervals)//2):
            #    self.__count_arrays[area][intervals[i]:intervals[i+1]] += 1

    def set_areas_value(self, areas, value):
        for area, intervals in areas.items():
            for i in range(len(intervals)//2):
                self.__count_arrays[area][intervals[i]:intervals[i+1]] = value

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
            self.__count_arrays[region_path][start:end] = value

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

    def find_valued_areas(self, value=False):
        # [1 1 1 2 2 2 0 0 0 1 1 1]
        # [- 0 0 1 0 0 1]
        areas = {}
        for node_id, count_array in self.__count_arrays.items():
            diffs = count_array[1:]-count_array[:-1]
            changes = np.where(diffs != 0)[0]
            changes = np.pad(changes, 1, "constant")
            changes[0] = -1
            changes += 1
            changes[-1] = self.graph.node_size(node_id)
            vals = count_array[changes[:-1]]
            start_idxs = np.where(vals == value)[0]
            end_idxs = start_idxs + 1
            starts = changes[start_idxs]
            ends = changes[end_idxs]
            areas[node_id] = [int(i) for i in chain.from_iterable(zip(starts, ends))]
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


class SmallIntervals(object):
    def __init__(self, interval_dict, end_interval_dict,
                 whole_intervals, graph, max_size):
        self.graph = graph
        self.interval_dict = interval_dict
        self.end_interval_dict = end_interval_dict
        self.whole_intervals = whole_intervals
        self.max_size = max_size

    def run(self):
        self.merge_intervals()
        self.filter_intervals()

    def merge_intervals(self):
        self.intervals = []
        for node_id, interval in self.end_interval_dict.items():
            intervals = self._get_intervals(node_id, [interval[0]])
            self.intervals.extend(
                [obg.Interval(i[0], i[-1], i[1:-1], graph=self.graph)
                 for i in intervals])

    def get_normal_intervals(self):
        self.small_areas = {}
        real_intervals = []
        for node_id, intervals in self.interval_dict.items():
            for i in range(len(intervals)//2):
                start = intervals[i]
                end = intervals[i+1]
                if start == 0 and all(
                        prev_node in self.end_interval_dict or prev_node in self.whole_intervals
                        for prev_node in
                        self.graph.reverse_adj_list[node_id]):
                    continue
                real_intervals.append(
                    Interval(start, end, [node_id], graph=self.graph))
        return real_intervals

    def filter_intervals(self):
        self.small_intervals = [interval for interval in self.intervals if interval.length()<=self.max_size]
        self.small_areas = {}
        for node_id, intervals in self.interval_dict.items():
            cur_list = []
            for i in range(len(intervals)//2):
                start = intervals[i]
                end = intervals[i+1]
                if start == 0 and all(
                        prev_node in self.end_interval_dict or prev_node in self.whole_intervals
                        for prev_node in
                        self.graph.reverse_adj_list[node_id]):
                    continue
                cur_list.extend([start, end])
            if cur_list:
                self.small_areas[node_id] = cur_list

    def filter_small_intervals(self):
        self.small_intervals = [interval for interval in self.intervals if interval.length()>=self.min_size]
        self.small_areas = {}
        for node_id, intervals in self.interval_dict.items():
            cur_list = []
            for i in range(len(intervals)//2):
                start = intervals[i]
                end = intervals[i+1]
                if start == 0 and all(
                        prev_node in self.end_interval_dict or prev_node in self.whole_intervals
                        for prev_node in
                        self.graph.reverse_adj_list[node_id]):
                    continue
                cur_list.extend([start, end])
            if cur_list:
                self.small_areas[node_id] = cur_list

    def _get_intervals(self, node_id, cur_interval):
        intervals = []
        my_interval = cur_interval + [node_id]
        for next_node in self.graph.adj_list[node_id]:
            if next_node in self.whole_intervals:
                intervals.extend(
                    self._get_intervals(next_node, my_interval))
                continue

            if next_node in self.interval_dict:
                if self.interval_dict[next_node] and self.interval_dict[next_node][0][0] == 0:
                    intervals.append(my_interval + [next_node, self.interval_dict[next_node][0][1]])
                    continue
            intervals.append(my_interval + [self.graph.node_size(node_id)])
        return intervals


class SparsePileup(Pileup):
    def create_data_struct(self):
        # valued intervals = [(idx,  value)]
        self.valued_intervals = defaultdict(list)

    def _add_areas(self, areas):
        for node_id, intervals in areas:
            old_intervals = self.valued_intervals
