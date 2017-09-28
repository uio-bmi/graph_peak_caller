from offsetbasedgraph.graphtraverser import GraphTraverser, update_areas
from offsetbasedgraph import Interval


def area_from_interval(interval, graph):
    areas = {}
    for i, region_path in enumerate(interval.region_paths):
        start = 0
        try:
            end = graph.node_size(region_path)
        except:
            print(interval)
            raise
        if region_path == interval.start_position.region_path_id:
            start = interval.start_position.offset
        if region_path == interval.end_position.region_path_id:
            end = interval.end_position.offset
        areas[region_path] = [start, end]
    return areas


def update_areas_fast(areas, new_areas):
    for area, startend in new_areas.items():
        if area not in areas:
            areas[area] = startend
            continue
        areas[area] = [min(startend[0], areas[area][0]),
                       max(startend[1], areas[area][1])]


class Shifter(object):
    def __init__(self, graph, d):
        self.graph = graph
        self.traverser = GraphTraverser(graph)
        self.d = d
        self.direction = self.d//abs(self.d) if self.d != 0 else 0

    def shift_interval(self, interval):
        start_positions = self.traverser.guided_shift(interval, self.d)
        length = interval.length()*self.direction
        areas = {}
        if interval.contains_position(start_positions[0]):
            if self.d < 0:
                trunc_interval = Interval(
                    interval.start_position, start_positions[0],
                    interval.region_paths[
                        0:interval.region_paths.index(
                            start_positions[0].region_path_id)])
            else:
                trunc_interval = Interval(
                    start_positions[0], interval.end_position,
                    interval.region_paths[
                        interval.region_paths.index(
                            start_positions[0].region_path_id):])

            areas = area_from_interval(trunc_interval, self.graph)
            start_positions = [interval.end_position] if self.d >= 0 else [interval.start_position]
            length = self.d

        for start_position in start_positions:
            new_areas = self.traverser.get_areas_from_point(
                start_position, length)

            update_areas(areas, new_areas)

        return areas

    def extend_interval(self, interval, local_direction=1):
        """Direction: -1 backward, +1 forward, 0 both"""
        assert local_direction in [-1, 0, 1]
        interval.graph = self.graph
        extension_length = self.d - interval.length()
        assert extension_length > 0
        areas = area_from_interval(interval, self.graph)
        if local_direction != -1:
            direction = interval.direction
            end_position = interval.end_position if direction == 1 else interval.start_position
            new_areas = self.traverser.get_areas_from_point(
                end_position, direction*extension_length)
            update_areas(areas, new_areas)
        if local_direction != 1:
            direction = interval.direction * -1
            end_position = interval.end_position if direction == 1 else interval.start_position
            new_areas = self.traverser.get_areas_from_point(
                end_position, direction*self.d)
            update_areas(areas, new_areas)

        return areas

    def get_areas_from_point(self, point, is_reverse, length):
        if is_reverse:
            length = (self.graph.node_size(point.region_path_id)-point.offset)+length
        else:
            length = point.offset + length
        visited = {}
        self.traverser.extend_from_block(
            point.region_path_id, length, visited)
        visited = {node_id: min(self.graph.node_size(node_id), l)
                   for node_id, l in visited.items()}
        if is_reverse:
            visited = {node_id: [self.graph.node_size(node_id)-l,
                                 self.graph.node_size(node_id)]
                       for node_id, l in visited.items()}
            visited[point.region_path_id][1] = point.offset
        else:
            visited = {node_id: [0, l] for node_id, l in visited.items()}
            visited[point.region_path_id][0] = point.offset
        return visited

    def extend_interval_fast(self, interval, local_direction=1):
        """Direction: +1 forward, 0 both"""
        assert local_direction in [-1, 0, 1]
        interval.graph = self.graph
        extension_length = self.d - interval.length()
        assert extension_length > 0, "Extension length is %d. D: %d, interval length: %d" % (extension_length, self.d, interval.length())
        areas = area_from_interval(interval, self.graph)
        direction = interval.direction
        end_position = interval.end_position if direction == 1 else interval.start_position
        new_areas = self.get_areas_from_point(
            end_position, direction == -1, extension_length)
        update_areas_fast(areas, new_areas)

        if local_direction == 0:
            direction = interval.direction * -1
            end_position = interval.end_position if direction == 1 else interval.start_position
            new_areas = self.get_areas_from_point(
                end_position, direction == -1, self.d)
            update_areas_fast(areas, new_areas)
        return areas
