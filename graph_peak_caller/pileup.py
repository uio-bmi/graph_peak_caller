import numpy as np


class Pileup(object):
    def __init__(self, graph, intervals):
        self.graph = graph
        self.intervals = intervals


    def __eq__(self, other):
        if False and self.graph != other.graph:
            return False
        for key, value in self.count_arrays.items():
            if key not in other.count_arrays:
                return False
            if not all(value == other.count_arrays[key]):
                return False

        return len(self.count_arrays) == len(other.count_arrays)

    @classmethod
    def from_bed_graph(cls, graph, file_name):
        file = open(file_name)
        pileup = cls(graph, [])
        pileup.create_count_arrays()

        graph_uses_int_ids = isinstance(list(graph.blocks.keys())[0], int)

        for line in file:
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
        self.count_arrays = {
            node_id: value*np.ones(block.length(), dtype="float")
            for node_id, block in self.graph.blocks.items()}

    def create_count_arrays(self):
        self.count_arrays = {node_id: np.zeros(block.length(), dtype="float")
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

    def add_area(self, block_id, start, end, value=1):
        count_array_block = self.count_arrays[block_id]
        count_array_block[start:end] += value

    def add_areas(self, areas):
        for area, intervals in areas.items():
            for i in range(len(intervals)//2):
                self.count_arrays[area][intervals[i]:intervals[i+1]] += 1

    def to_bed_graph(self, filename):
        f = open(filename, "w")
        for node_id, count_array in self.count_arrays.items():
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
        return filename

    def summary(self):
        return sum(array.sum() for array in self.count_arrays.values())


class SparsePileup(Pileup):
    def create_data_struct(self):
        # valued intervals = [(idx,  value)]
        self.valued_intervals = defaultdict(list)

    def _add_areas(self, areas):
        for node_id, intervals in areas:
            old_intervals = self.valued_intervals
