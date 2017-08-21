import numpy as np
from shifter import Shifter


class Pileup(object):
    def __init__(self, graph, intervals, shift=0):
        self.graph = graph
        self.intervals = intervals
        self.shifter = Shifter(graph, intervals, shift)
        self.shift = shift

    def __eq__(self, other):
        if False and self.graph != other.graph:
            return False
        for key, value in self.count_arrays.items():
            if key not in other.count_arrays:
                return False
            if not all(value == other.count_arrays[key]):
                return False

        return len(self.count_arrays) == len(other.count_arrays)

    def create(self):
        print("Creating")
        self.create_count_arrays()
        for interval in self.intervals:
            if self.shift != 0:
                self.add_extended_interval(interval)
            else:
                self.add_interval(interval)

    def init_value(self, value):
        self.count_arrays = {
            node_id: value*np.ones(block.length(), dtype="int32")
            for node_id, block in self.graph.blocks.items()}

    def create_count_arrays(self):
        print("Init count_arrays")
        self.count_arrays = {node_id: np.zeros(block.length(), dtype="int32")
                             for node_id, block in self.graph.blocks.items()}

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
            self.count_arrays[region_path][start:end] += 1

    def add_shifted_interval(self, interval):
        areas = self.shifter.shift_interval(interval)
        for area, intervals in areas.items():
            for i in range(len(intervals)//2):
                self.count_arrays[area][intervals[i]:intervals[i+1]] += 1

    def add_extended_interval(self, interval):
        areas = self.shifter.extend_interval(interval)
        for area, intervals in areas.items():
            for i in range(len(intervals)//2):
                self.count_arrays[area][intervals[i]:intervals[i+1]] += 1

    def to_bed_graph(self, filename):
        f = open(filename, "w")
        for node_id, count_array in self.count_arrays.items():
            start = None
            cur_val = count_array[0]
            for i, new_val in enumerate(count_array):
                if new_val != cur_val:
                    if start is not None:
                        interval = (node_id, start, i, cur_val)
                        f.write("%s\t%s\t%s\t%s\n" % interval)
                    start = i
                    cur_val = new_val
            interval = (node_id, start, len(count_array), cur_val)
            f.write("%s\t%s\t%s\t%s\n" % interval)

    def summary(self):
        return sum(array.sum() for array in self.count_arrays.values())


class SparsePileup(Pileup):
    def create_data_struct(self):
        # valued intervals = [(idx,  value)]
        self.valued_intervals = defaultdict(list)

    def _add_areas(self, areas):
        for node_id, intervals in areas:
            old_intervals = self.valued_intervals
            for node_id
