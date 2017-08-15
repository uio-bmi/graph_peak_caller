import numpy as np


class Pileup(object):
    def __init__(self, graph, intervals):
        self.graph = graph
        self.intervals = intervals

    def create(self):
        self.create_count_arrays()
        for interval in self.intervals:
            self.add_interval(interval)

    def create_count_arrays(self):
        self.count_arrays = {node_id: np.zeros(block.length(), dtype="int32")
                             for node_id, block in self.graph.blocks.items()}

    def add_interval(self, interval):
        if not all(region_path in self.graph.blocks for region_path in interval.region_paths):
            return
        print("Adding")
        for i, region_path in enumerate(interval.region_paths):
            start = 0
            end = self.graph.blocks[region_path].length()
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths)-1:
                end = interval.end_position.offset
            print(region_path, start, end, self.graph.blocks[region_path].length())
            self.count_arrays[region_path][start:end] += 1

    def summary(self):
        return sum(array.sum() for array in self.count_arrays.values())
