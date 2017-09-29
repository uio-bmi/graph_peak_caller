import logging

from offsetbasedgraph.graphtraverser import GraphTraverser
from offsetbasedgraph.interval import Position
import offsetbasedgraph as obg
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
        f = all if include_partial_stubs else any
        intervals = []
        for node_id, startends in self.areas.items():
            for i in range(len(startends)//2):
                start = startends[2*i]
                end = startends[2*i+1]
                if start == 0:
                    if f(prev_id in self.areas and self.areas[prev_id] and
                         self.areas[prev_id][-1] == self.graph.node_size(prev_id)
                         for prev_id in self.graph.reverse_adj_list[node_id]):
                        continue
                if end == self.graph.node_size(node_id):
                    intervals.extend(
                        self._get_intervals(
                            node_id, [start],
                            self.areas, include_partial_stubs))
                else:
                    intervals.append(
                        obg.DirectedInterval(int(start), int(end),
                                             [node_id], graph=self.graph))
        return intervals


class OldExtender(object):
    def __init__(self, graph, d):
        self.graph = graph
        self.pos_traverser = GraphTraverser(graph)
        self.neg_traverser = GraphTraverser(graph, -1)
        logging.debug(self.pos_traverser)
        logging.debug(self.neg_traverser)
        self.d = d

    def get_areas_from_point(self, point, length, traverser):
        length = point.offset + length
        visited = {}
        traverser.extend_from_block(
            point.region_path_id, length, visited)
        visited = {node_id: min(self.graph.node_size(node_id), l)
                   for node_id, l in visited.items()}
        visited = {node_id: [0, l] for node_id, l in visited.items()}
        visited[point.region_path_id][0] = point.offset
        return Areas(self.graph, visited)

    def get_areas_from_point_new(self, point, length, traverser):
        inverse_start = self.graph.node_size(point.region_path_id)
        remaining = length-inverse_start
        if remaining < 0:
            return Areas(self.graph, {point.node_id:
                                      [point.offset, point.offset+length]})
        visited = {-point.region_path_id: [0, inverse_start]}
        for next_node in traverser.adj_list[point.region_path_id]:
            traverser.extend_from_block(next_node, remaining, visited)
        visited = {node_id: min(self.graph.node_size(node_id), l)
                   for node_id, l in visited.items()}
        visited = {node_id: [0, l] for node_id, l in visited.items()}
        # visited[point.region_path_id][0] = point.offset
        return Areas(self.graph, visited)

    def extend_interval(self, interval, local_direction=1):
        """Direction: +1 forward, 0 both"""
        assert local_direction in [-1, 0, 1]
        interval.graph = self.graph
        extension_length = self.d - interval.length()
        areas = Areas.from_interval(interval, self.graph)
        end_position = interval.end_position
        if interval.can_be_on_positive_strand():
            new_areas = self.get_areas_from_point(
                end_position, extension_length, self.pos_traverser)
            areas.update(new_areas)
        if interval.can_be_on_negative_strand():
            logging.debug("NEGATIVE")
            new_areas = self.get_areas_from_point(
                end_position, extension_length, self.neg_traverser)
            logging.info("New areas: %s" % new_areas)
            areas.update(new_areas)

        if local_direction == 0:
            end_position = interval.start_position
            region_path = end_position.region_path_id * (-1)
            offset = self.graph.node_size(region_path) - end_position.offset
            end_position = Position(region_path, offset)

            if interval.can_be_on_positive_strand():
                new_areas = self.get_areas_from_point(
                    end_position, self.d, self.neg_traverser)
                areas.update(new_areas)
            if interval.can_be_on_negative_strand():
                new_areas = self.get_areas_from_point(
                    end_position, self.d, self.pos_traverser)
                areas.update(new_areas)
        areas.reverse_reversals()
        logging.warning("Extending")
        logging.warning(interval)
        logging.warning(areas)
        return areas


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
        self.graph = graph
        self.pos_traverser = GraphTraverser(graph, +1)
        self.neg_traverser = GraphTraverser(graph, -1)

    def get_areas_from_node(self, region_path, length, traverser):
        logging.debug("########")
        logging.debug(region_path)
        logging.debug(length)
        visited = {}
        for next_node in traverser.adj_list[region_path]:
            traverser.extend_from_block(next_node, length, visited)
        visited = {node_id: min(self.graph.node_size(node_id), l)
                   for node_id, l in visited.items()}
        logging.debug(visited)
        self.area_builder.update(
             {node_id: [0, l] for node_id, l in visited.items()})

    def extend_interval(self, interval, direction=1):
        self.area_builder = AreasBuilder(self.graph)
        logging.warning(interval.length())
        pos_length = self.length-interval.length()
        neg_length = self.length if direction == 0 else 0
        pos_remain, neg_remain = self.area_builder.filled_interval(
            interval, pos_length, neg_length)
        logging.debug(pos_remain)
        logging.debug(neg_remain)
        logging.debug(self.area_builder.areas)
        is_pos = interval.can_be_on_positive_strand()
        is_neg = interval.can_be_on_negative_strand()
        logging.debug(is_pos)
        logging.debug(is_neg)
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

        self.area_builder.reverse_reversals()
        return Areas(self.graph, self.area_builder.areas)
