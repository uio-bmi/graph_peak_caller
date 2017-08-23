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
        return "\n".join("%s: %s" % (node_id, sum(array)) for node_id, array in self.__count_arrays.items())

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

    def threshold(self, cutoff):
        for node_id, count_array in self.__count_arrays.items():
            self.__count_arrays[node_id] = count_array > cutoff

    def clear(self):
        self.__count_arrays = None

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
            end = self.graph.blocks[region_path].length()
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
            for i in range(len(intervals)//2):
                self.__count_arrays[area][intervals[i]:intervals[i+1]] += 1

    def set_areas_value(self, areas, value):
        for area, intervals in areas.items():
            for i in range(len(intervals)//2):
                self.__count_arrays[area][intervals[i]:intervals[i+1]] = value

    def set_interval_value(self, interval, value):
        assert all(region_path in self.graph.blocks for
                   region_path in interval.region_paths)
        for i, region_path in enumerate(interval.region_paths):
            start = 0
            end = self.graph.blocks[region_path].length()
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths)-1:
                end = interval.end_position.offset
            self.__count_arrays[region_path][start:end] = value

    def to_bed_graph(self, filename):
        f = open(filename, "w")
        for node_id, count_array in self.__count_arrays.items():
            start = 0
            cur_val = None
            for i, new_val in enumerate(count_array):
                if new_val != cur_val:
                    if cur_val is not None:
                        interval = (node_id, start, i, cur_val)
                        f.write("%s\t%s\t%s\t%s\n" % interval)
                    start = i
                    cur_val = new_val
            interval = (node_id, start, len(count_array), cur_val)
            f.write("%s\t%s\t%s\t%s\n" % interval)
        f.close()
        self.filename = filename
        self.is_written = True
        return filename

    def scale(self, scale):
        for count_array in self.__count_arrays.values():
            count_array *= scale

    def summary(self):
        return sum(array.sum() for array in self.__count_arrays.values())

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
        interval_dict, end_interval_dict, whole_intervals = self.find_valued_intevals(0)
        small_intervals = SmallIntervals(interval_dict, end_interval_dict, whole_intervals, self.graph, max_size)
        small_intervals.run()
        self.set_areas_value(small_intervals.small_areas, 1)
        [self.set_interval_value(interval, 1) for interval in small_intervals.small_intervals]


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

    def filter_intervals(self):
        self.small_intervals = [interval for interval in self.intervals if interval.length()<=self.max_size]
        self.small_areas = {}
        for node_id, intervals in self.interval_dict.items():
            cur_list = []
            for i in range(len(intervals)//2):
                start = intervals[i]
                end = intervlas[i+1]
                if start == 0 and all(prev_node in self.end_interval_dict for prev_node in
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
            intervals.append(my_interval + [self.graph.blocks[node_id].length()])
        return intervals


class SparsePileup(Pileup):
    def create_data_struct(self):
        # valued intervals = [(idx,  value)]
        self.valued_intervals = defaultdict(list)

    def _add_areas(self, areas):
        for node_id, intervals in areas:
            old_intervals = self.valued_intervals
