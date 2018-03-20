import numpy as np


class DensePileup:
    def __init__(self, graph, values):
        self._values = values
        self._graph = graph

    def values_in_range(self, node, start, end):
        index = node - self._graph.min_node
        array_start = self._graph.node_indexes[index] + start
        array_end = self._graph.node_indexes[index] + end
        assert array_end > array_start
        assert array_end <= len(self._values), "Array end %d > len of values %s" % (array_end, len(self._values))
        out = self._values[array_start:array_end]
        assert len(out) > 0
        return out

    def values(self, node_id):
        index = node_id - self._graph.min_node
        array_start = self._graph.node_indexes[index]
        array_end = self._graph.node_indexes[index+1]
        return self._values[array_start:array_end]

    def get_interval_values(self, interval):
        interval_length = interval.length()
        assert interval_length > 0, "Trying to get value of interval with negative length, %s" % interval
        values = np.zeros(interval_length)
        offset = 0

        is_reverse = False
        if interval.region_paths[0] < 0:
            assert np.all(np.array(interval.region_paths) < 0), " First region path negative, but not rest. Interval: %s" % interval
            is_reverse = True
            interval = interval.get_reverse()

        for i, rp in enumerate(interval.region_paths):
            assert rp > 0, "Currently only implemented for forward directed intervals"
            start = 0
            end = self._graph.node_size(rp)
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths) - 1:
                end = interval.end_position.offset

            assert end > start, "Start >= end for interval %s" % interval
            values_in_rp = self.values_in_range(rp, start, end)
            values[offset:offset + (end - start)] = values_in_rp

            offset += end-start


        if is_reverse:
            out = values[::-1]
        else:
            out = values

        assert np.sum(np.isnan(values)) == 0, "%s contains nan" % (values)
        return values

    def __str__(self):
        out = "Minimum DensePileup. Values: %s. Sum: %.3f" % (self._values, np.sum(self._values))
        return out

    def __repr__(self):
        return out
