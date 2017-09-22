from offsetbasedgraph.graphtraverser import GraphTraverser
from offsetbasedgraph.interval import Position


class Areas(object):
    def __init__(self, graph, areas=None):
        if areas is None:
            self.areas = {}
        else:
            self.areas = areas

    @classmethod
    def from_interval(cls, interval, graph):
        areas = {}
        for region_path in interval.region_paths:
            start = 0
            end = graph.node_size(region_path)
            if region_path == interval.start_position.region_path_id:
                start = interval.start_position.offset
            if region_path == interval.end_position.region_path_id:
                end = interval.end_position.offset
            areas[region_path] = [start, end]

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

    def reverse_reversals(self):
        for node_id, startend in self.areas.items():
            if node_id >= 0:
                continue
            l = self.graph.node_size(node_id)
            pos_coords = [l-pos for pos in reversed(startend)]
            if -node_id not in self.areas:
                self.areas[-node_id] = pos_coords
            forward_areas = self.areas[-node_id]
            if forward_areas[-1] < pos_coords[0]:
                forward_areas.extend(pos_coords)
            else:
                self.areas[-node_id] = [forward_areas[0], pos_coords[-1]]
        self.areas = {key: value for key, value in self.areas.items()
                      if key >= 0}

    def get_starts(self, node_id):
        return [idx for i, idx in enumerate(self.areas[node_id])
                if i % 2 == 0]

    def get_ends(self, node_id):
        return [idx for i, idx in enumerate(self.areas[node_id])
                if i % 2 == 1]


class Extender(object):
    def __init__(self, graph, d):
        self.graph = graph
        self.pos_traverser = GraphTraverser(graph)
        self.neg_traverser = GraphTraverser(graph, -1)
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

    def extend_interval_fast(self, interval, local_direction=1):
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
            new_areas = self.get_areas_from_point(
                end_position, extension_length, self.neg_traverser)
            areas.update(new_areas)

        if local_direction == 0:
            end_position = interval.start_position
            region_path = end_position.region_path * (-1)
            offset = self.graph.node_size(region_path) - end_position.offset
            end_position = Position(region_path, offset)

            if interval.can_be_on_positive_strand():
                new_areas = self.get_areas_from_point(
                    end_position, self.d, self.neg_traverser)
                areas.update(areas, new_areas)
            if interval.can_be_on_negative_strand():
                new_areas = self.get_areas_from_point(
                    end_position, self.d, self.pos_traverser)
                areas.update(areas, new_areas)
        areas.reverse_reversals()
        return areas
