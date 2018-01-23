import logging
from collections import defaultdict
from offsetbasedgraph.graphtraverser import GraphTraverser
from offsetbasedgraph.interval import Position
from .areas import BinaryContinousAreas
import offsetbasedgraph as obg
import numpy as np
import itertools
# logging.basicConfig(level=logging.DEBUG)


class Areas(object):
    def __init__(self, graph, areas=None):
        self.graph = graph
        if areas is None:
            self.areas = {}
        else:
            self.areas = areas

    def __eq__(self, other):
        if not len(self.areas) == len(other.areas):
            return False

        if not all(node_id in other.areas for node_id in self.areas):
            return False

        return all(startend == other.areas[node_id] for node_id, startend
                   in self.areas.items())

    def __str__(self):
        return str(self.areas)

    def __repr__(self):
        return repr(self.areas)

    def join(self, areas):
        pass

    @classmethod
    def from_interval(cls, interval, graph):
        if len(interval.region_paths) == 1:
            return cls(graph, {interval.region_paths[0]: [
                interval.start_position.offset,
                interval.end_position.offset]})
        areas = {}
        start_rp = interval.region_paths[0]
        start_len = graph.node_size(start_rp)-interval.start_position.offset
        areas[-start_rp] = [0, start_len]
        areas[interval.end_position.region_path_id] = [0, interval.end_position.offset]

        for region_path in interval.region_paths[1:-1]:
            areas[region_path] = [0, graph.node_size(region_path)]

        return cls(graph, areas)

    def update(self, other):
        """TODO: Hanlde cyclic graphs"""
        for node_id, startend in other.areas.items():
            if node_id not in self.areas:
                self.areas[node_id] = startend
                continue
            self.areas[node_id] = [
                min(startend[0], self.areas[node_id][0]),
                max(startend[1], self.areas[node_id][1])]

    def _add_nontrivial_areas_for_node(self, node_id, starts_and_ends):
        raise NotImplementedError("Nontrivial addition not supported")

    def add_areas_for_node(self, node_id, starts_and_ends):
        if node_id not in self.areas:
            self.areas[node_id] = starts_and_ends
        else:
            current = self.areas[node_id]
            if len(current) > 2 or len(starts_and_ends) > 2:
                return self._add_nontrivial_areas_for_node(
                    node_id, starts_and_ends)

            if current[1] < starts_and_ends[0]:
                new = np.append(current, starts_and_ends)
            elif current[0] > starts_and_ends[1]:
                new = np.append(starts_and_ends, current)
            else:
                # todo
                print("Error")
                print(current)
                print(starts_and_ends)
                raise NotImplementedError(
                    "Case where new start equals old end is not implemented")

            assert len(new) == 4
            self.areas[node_id] = new

    def robust_update(self, other):
        for node_id, startend in other.areas.items():
            if node_id not in self.areas:
                self.areas[node_id] = startend
                continue

            self.areas[node_id] = [
                min(startend[0], self.areas[node_id][0]),
                max(startend[1], self.areas[node_id][1])]

    def reverse_reversals(self):
        neg_node_ids = [node_id for node_id in self.areas.keys()
                        if node_id < 0]
        for node_id in neg_node_ids:
            
            startend = self.areas[node_id]
            l = self.graph.node_size(node_id)
            pos_coords = [l-pos for pos in reversed(startend)]
            del self.areas[node_id]
            if -node_id not in self.areas:
                self.areas[-node_id] = pos_coords
                continue
            forward_areas = self.areas[-node_id]

            if forward_areas[-1] < pos_coords[0]:
                forward_areas.extend(pos_coords)
            else:
                all_coords = forward_areas + pos_coords
                self.areas[-node_id] = [min(all_coords), max(all_coords)]

    def get_starts(self, node_id):
        return [idx for i, idx in enumerate(self.areas[node_id])
                if i % 2 == 0]

    def get_ends(self, node_id):
        return [idx for i, idx in enumerate(self.areas[node_id])
                if i % 2 == 1]

    def to_intervals(self, include_partial_stubs):
        raise NotImplementedError()

    def to_simple_intervals(self):
        intervals = []
        for node, area in self.areas.items():
            for i in range(0, len(area) // 2):
                start = int(area[i*2])
                end = int(area[i*2+1])
                interval = obg.Interval(start, end, [node], self.graph)
                intervals.append(interval)

        return intervals

    def _is_start_or_end_position(self, node, position_offset):
        if self._is_end_position(node, position_offset) \
                or self._is_start_position(node, position_offset):
            return True

        return False

    def _is_end_position(self, node, offset):
        if offset < self.graph.node_size(node):
            return True
        graph = self.graph
        for in_node in itertools.chain(graph.adj_list[node], graph.reverse_adj_list[node]):
            if in_node in self.areas and len(self.areas[in_node]) > 0:
                if self.areas[in_node][0] == 0:
                    return False
            elif -in_node in self.areas and len(self.areas[-in_node]) > 0:
                if self.areas[-in_node][-1] == graph.node_size(in_node):
                    return False

        return True

    def _is_start_position(self, node, offset):
        if offset > 0:
            return True

        graph = self.graph
        for in_node in itertools.chain(graph.adj_list[-node], graph.reverse_adj_list[-node]):
            if in_node in self.areas and len(self.areas[in_node]) > 0:
                if self.areas[in_node][0] == 0:
                    return False

            elif -in_node in self.areas and len(self.areas[-in_node]) > 0:
                if self.areas[-in_node][-1] == graph.node_size(in_node):
                    return False

        return True

    def get_start_and_end_positions(self):
        positions = []
        for node, starts_and_ends in self.areas.items():
            n_starts_ends = len(starts_and_ends)
            for start in starts_and_ends[range(0, n_starts_ends, 2)]:
                if self._is_start_position(node, start):
                    positions.append(obg.Position(node, start))

            for end in starts_and_ends[range(1, n_starts_ends, 2)]:
                if self._is_end_position(node, end):
                    positions.append(obg.Position(node, end-1))
                    # Correct to non-inclusive end, since this may be interpreted as a start

        return positions

    def get_all_included_nodes(self):
        return self.areas.keys()

    def to_file_line(self):
        start_ends = ','.join([str(position) for position in
                               self.get_start_and_end_positions()])
        nodes = ','.join([str(node) for node in self.get_all_included_nodes()])
        return "%s\t%s\n" % (start_ends, nodes)

    @classmethod
    def from_file_line(cls, line):
        pass


class AreasBuilder(object):
    def __init__(self, graph):
        self.graph = graph
        self.areas = {}

    def update(self, new_areas):
        """NB!: Only works if all intervals start at 0"""
        for rp, startend in new_areas.items():
            if rp not in self.areas:
                self.areas[rp] = startend
                continue
            self.areas[rp][-1] = max(startend[-1], self.areas[rp][-1])

    def reverse_reversals(self):
        neg_rps = [rp for rp in self.areas if rp < 0]
        for rp in neg_rps:
            node_size = self.graph.node_size(rp)
            pos_rp = -rp
            if pos_rp not in self.areas:
                self.areas[pos_rp] = [
                    node_size-i for i in self.areas[rp][::-1]]
            else:
                mid = node_size-self.areas[rp][-1]
                if self.areas[pos_rp][-1] > mid:
                    self.areas[pos_rp] = [0, node_size]
                else:
                    self.areas[pos_rp] += [mid, node_size]
            del self.areas[rp]

    def filled_interval(self, interval, fill_pos=0, fill_neg=0):
        start = max(interval.start_position.offset-fill_neg, 0)
        end_size = self.graph.node_size(interval.end_position.region_path_id)
        end = min(interval.end_position.offset+fill_pos,
                  end_size)
        neg_remain = max(0, fill_neg-interval.start_position.offset)
        pos_remain = max(0, fill_pos-(end_size-interval.end_position.offset))

        if len(interval.region_paths) == 1:
            # Make sure it touches edge if possible
            if start == 0:
                self.areas[interval.region_paths[0]] = [start, end]
            else:
                self.areas[-interval.region_paths[0]] = [end_size-end, end_size-start]
            return pos_remain, neg_remain
        start_rp = interval.region_paths[0]
        start_len = self.graph.node_size(start_rp) - start
        self.areas[-start_rp] = [0, start_len]
        self.areas[interval.end_position.region_path_id] = [0, end]

        for region_path in interval.region_paths[1:-1]:
            self.areas[region_path] = [0, self.graph.node_size(region_path)]

        return pos_remain, neg_remain


class Extender(object):
    def __init__(self, graph, length):
        self.length = length
        assert isinstance(length, int)
        self.graph = graph
        self.pos_traverser = GraphTraverser(graph, +1)
        self.neg_traverser = GraphTraverser(graph, -1)

    def get_areas_from_node(self, region_path, length, traverser):
        visited = defaultdict(int)
        for next_node in traverser.adj_list[region_path]:
            traverser.extend_from_block(next_node, length, visited)

        for node_id, l in visited.items():
            if l == 0:
                continue
            if l >= self.graph.node_size(node_id):
                self.area_builder.add_full(node_id)
            else:
                self.area_builder.add_start(node_id, l)

    def extend_interval(self, interval, direction=1):
        self.area_builder = BinaryContinousAreas(self.graph)
        pos_length = self.length-interval.length()
        neg_length = self.length if direction == 0 else 0
        pos_remain, neg_remain = self.area_builder.filled_interval(
            interval, pos_length, neg_length)
        is_pos = interval.can_be_on_positive_strand()
        is_neg = interval.can_be_on_negative_strand()
        #print(" Extending interval %s. Can be on pos/neg: %d/%d" % (interval, is_pos, is_neg))
        if pos_remain:
            if is_pos:
                self.get_areas_from_node(
                    interval.end_position.region_path_id,
                    pos_remain, self.pos_traverser)
            if is_neg:
                self.get_areas_from_node(
                    interval.end_position.region_path_id,
                    pos_remain, self.neg_traverser)
        if neg_remain:
            if is_pos:
                self.get_areas_from_node(
                    -interval.start_position.region_path_id,
                    neg_remain, self.neg_traverser)
            if is_neg:
                self.get_areas_from_node(
                    -interval.start_position.region_path_id,
                    neg_remain, self.pos_traverser)

        self.area_builder.sanitize()
        # self.area_builder.reverse_reversals()
        return self.area_builder
# Areas(self.graph, self.area_builder.areas)
