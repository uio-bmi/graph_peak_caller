import numpy as np


class DensePileup:
    def __init__(self, graph, values):
        self._values = values
        self._graph = graph

    def values_in_range(self, node, start, end):
        index = node - self._graph.min_node
        array_start = self._graph.node_indexes[index] + start
        array_end = self._graph.node_indexes[index] + end
        return self._values[array_start:array_end]

    def values(self, node_id):
        index = node_id - self._graph.min_node
        array_start = self._graph.node_indexes[index]
        array_end = self._graph.node_indexes[index+1]
        return self._values[array_start:array_end]

    def get_interval_values(self, interval):
        values = np.zeros(interval.length())
        offset = 0
        if all([rp < 0 for rp in interval.region_paths]):
            # Reverse before finding
            # find_reversed = True
            interval = interval.get_reverse()

        for i, rp in enumerate(interval.region_paths):
            assert rp > 0, "Currently only implemented for forward directed intervals"
            start = 0
            end = self._graph.node_size(rp)
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths) - 1:
                end = interval.end_position.offset

            values_in_rp = self.values_in_range(rp, start, end)
            values[offset:offset + (end - start)] = values_in_rp

            offset += end-start

        return values
