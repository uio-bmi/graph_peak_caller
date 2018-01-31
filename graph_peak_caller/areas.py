from itertools import chain
from collections import defaultdict
import numpy as np
import offsetbasedgraph as obg
import itertools
from .ioclass import CollectionIO

class Areas(object):
    pass


class BinaryContinousAreas(Areas):
    def __init__(self, graph):
        self.graph = graph
        self.full_areas = {}
        self.starts = defaultdict(int)
        self.internal_intervals = {}

    def __str__(self):
        full_str = "Full: %s" % ",".join(
            str(node_id) for node_id in self.full_areas)
        start_str = "Start: %s" % ", ".join(
            "(%s,%s)" % (node_id, idx) for node_id, idx in self.starts.items())
        internal_str = "Internsals: %s" % self.internal_intervals
        return "\n".join((full_str, start_str, internal_str))

    __repr__ = __str__

    def __eq__(self, other):
        if self.full_areas != other.full_areas:
            return False
        if self.starts != other.starts:
            return False
        return self.internal_intervals == other.internal_intervals

    def add_full(self, node_id):
        self.full_areas[abs(node_id)] = 1

    def add_start(self, node_id, idx):
        assert idx > 0
        self.starts[node_id] = max(idx, self.starts[node_id])

    def _is_end_position(self, node, offset):
        if offset < self.graph.node_size(node):
            return True
        graph = self.graph
        for in_node in itertools.chain(graph.adj_list[node], graph.reverse_adj_list[node]):
            if in_node in self.full_areas or in_node in self.starts or -in_node in self.full_areas:
                return False

        return True

    def _is_start_position(self, node, offset):
        if offset > 0:
            return True

        graph = self.graph
        for in_node in itertools.chain(graph.adj_list[-node], graph.reverse_adj_list[-node]):
            if in_node in self.full_areas or in_node in self.starts or -in_node in self.full_areas:
                return False

        return True

    def merge_with_other(self, other):
        for full in other.full_areas.keys():
            self.add_full(full)
        for node, start in other.starts.items():
            self.add_start(node, start)
        for node, startend in other.internal_intervals.items():
            self.add_internal(node, startend[0], startend[1])

    def add_internal(self, node_id, start, end):
        assert start != end
        self.internal_intervals[node_id] = [start, end]

    def add(self, node_id, start, end):
        assert start != end
        node_size = self.graph.node_size(node_id)
        if start == 0 and end == node_size:
            self.add_full(node_id)
        elif start == 0:
            self.add_start(node_id, end)
        elif end == node_size:
            self.add_start(-node_id, node_size-start)
        else:
            if node_id < 0:
                start, end = node_size-end, node_size-start
                node_id = -node_id
            self.add_internal(node_id, start, end)

    def join_starts_and_ends(self):
        negative_node_ids = [node_id for node_id in self.starts if
                             node_id < 0]
        for node_id in negative_node_ids:
            if -node_id not in self.starts:
                continue
            if self.starts[node_id] + self.starts[-node_id] < self.graph.node_size(-node_id):
                continue
            self.full_areas[-node_id] = 1
            del self.starts[node_id]
            del self.starts[-node_id]

    def sanitize(self):
        for node_id in self.full_areas:
            if node_id in self.starts:
                del self.starts[node_id]
            if -node_id in self.starts:
                del self.starts[-node_id]
        self.join_starts_and_ends()

    def _add_abs_internal(self, node_id, start, end):
        if node_id < 0:
            node_size = self.graph.node_size(node_id)
            rev_start = node_size-end
            rev_end = node_size-start
            self.add_internal(abs(node_id), rev_start, rev_end)
        else:
            self.add_internal(node_id, start, end)

    def filled_interval(self, interval, fill_pos=0, fill_neg=0):
        start = max(interval.start_position.offset-fill_neg, 0)
        end_size = self.graph.node_size(interval.end_position.region_path_id)
        end = min(interval.end_position.offset+fill_pos,
                  end_size)
        neg_remain = max(0, fill_neg-interval.start_position.offset)
        pos_remain = max(0, fill_pos-(end_size-interval.end_position.offset))

        if len(interval.region_paths) == 1:
            # Make sure it touches edge if possible
            rp = interval.region_paths[0]
            if start == 0 and end == end_size:
                self.add_full(abs(rp))
            elif start == 0:
                self.add_start(rp, end)
            elif end == end_size:
                self.add_start(-rp, end-start)
            else:
                self._add_abs_internal(rp, start, end)
            return pos_remain, neg_remain
        start_rp = interval.region_paths[0]
        start_len = self.graph.node_size(start_rp) - start
        self.add_start(-start_rp, start_len)
        self.add_start(interval.end_position.region_path_id, end)
        for region_path in interval.region_paths[1:-1]:
            self.add_full(abs(region_path))
        return pos_remain, neg_remain

    def get_start_positions(self):
        start_positions = [
            obg.Position(-node_id, self.graph.node_size(node_id)-offset)
            for node_id, offset in self.starts.items()]
        node_ids = list(self.full_areas.keys()) + [
            -node_id for node_id in self.full_areas]
        previous_nodes_list = [
            self.graph.reverse_adj_list[-node_id]
            for node_id in node_ids]

        def filter_my_nodes(node_list):
            return [node_id for node_id in node_list if
                    abs(node_id) in self.full_areas or node_id
                    in self.starts]
        previous_nodes_list = [filter_my_nodes(prev_nodes) for
                               prev_nodes in previous_nodes_list]
        full_starts = [obg.Position(node_id, 0) for node_id, prev_nodes in
                       zip(node_ids, previous_nodes_list)
                       if not prev_nodes]

        return start_positions + full_starts

    def get_node_ids(self):
        return chain(self.full_areas.keys(), self.starts.keys(),
                     self.internal_intervals.keys())

    @classmethod
    def from_old_areas(cls, areas):
        binary_connected_areas = cls(areas.graph)
        for node_id in areas.areas:
            starts = areas.get_starts(node_id)
            ends = areas.get_ends(node_id)
            for start, end in zip(starts, ends):
                binary_connected_areas.add(node_id, start, end)

        return binary_connected_areas

    def to_file_line(self):
        fulls = ",".join(str(node_id) for node_id in self.full_areas)
        starts = ",".join("%s:%s" % (node_id, offset)
                          for node_id, offset in self.starts.items())
        internals = ",".join("%s:%s-%s" % (node_id, v[0], v[1])
                             for node_id, v in self.internal_intervals.items())
        return "(%s)\t(%s)\t(%s)\n" % (fulls, starts, internals)

    @classmethod
    def from_file_line(cls, line, graph):
        obj = cls(graph)
        fulls, starts, internals = [v[1:-1] for v in line.split()]
        if fulls:
            obj.full_areas = {int(full): 1 for full in fulls.split(",")}
        if starts:
            obj.starts = {int(p.split(":")[0]): int(p.split(":")[1]) for
                          p in starts.split(",")}
        if internals:
            k, v = internals.split(":")
            obj.internal_intervals = {int(k): [int(v.split("-")[0]),
                                               int(v.split("-")[1])]}
        return obj


class BCACollection(CollectionIO):
    _obj_type = BinaryContinousAreas


class ValuedAreas(Areas):
    def __init__(self, graph):
        self.graph = graph
        self.full_areas = defaultdict(int)
        self.starts = defaultdict(list)
        self.internal_starts = defaultdict(list)
        self.internal_ends = defaultdict(list)

    def has_anything_on_node(self, node_id):
        if node_id in self.full_areas:
            return True
        if node_id in self.starts:
            return True
        if node_id in self.internal_ends:
            return True
        if node_id in self.internal_starts:
            return True
        return False

    def add_binary_areas(self, areas, touched_nodes=None):
        for node_id in areas.full_areas:
            self.full_areas[node_id] += 1
            touched_nodes.add(abs(node_id))
        for node_id, internal_intervals in areas.internal_intervals.items():
            self.internal_starts[node_id].append(internal_intervals[0])
            self.internal_ends[node_id].append(internal_intervals[1])
            touched_nodes.add(abs(node_id))
        for node_id, start in areas.starts.items():
            self.starts[node_id].append(start)
            touched_nodes.add(abs(node_id))

    def get_starts_array(self, node_id, node_size=None):
        if node_size is None:
            node_size = self.graph.node_size(node_id)
        """
        length_starts_full = len(self.starts[node_id])+self.full_areas[node_id]
        length_internal = len(self.internal_starts[node_id])
        length_starts = len(self.starts[-node_id])
        out = np.zeros(length_starts_full + length_internal + length_starts)
        b = length_starts_full+length_internal
        out[length_starts_full:b] = self.internal_starts[node_id]
        out[b:] = self.starts[-node_id]
        out[b:] *= -1
        out[b:] += node_size
        return out
        """

        starts = [0]*(len(self.starts[node_id])+self.full_areas[node_id])
        starts.extend(self.internal_starts[node_id])
        starts.extend([node_size-start for start in self.starts[-node_id]])
        return np.array(starts, dtype="int")

    def get_ends_array(self, node_id, node_size=None):
        if node_size is None:
            node_size = self.graph.node_size(node_id)

        ends = self.starts[node_id][:]
        ends.extend(self.internal_ends[node_id])
        ends.extend([node_size]*(
            len(self.starts[-node_id])+self.full_areas[node_id]))
        return np.array(ends, dtype="int")
