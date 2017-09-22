from offsetbasedgraph.graphtraverser import GraphTraverser
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


class Extender(object):
    def __init__(self, graph, d):
        self.graph = graph
        self.pos_traverser = GraphTraverser(graph)
        self.neg_traverser = GraphTraverser(graph, -1)
        self.d = d
        self.direction = self.d//abs(self.d) if self.d != 0 else 0

    def get_areas_from_point(self, point, length):
        length = point.offset + length
        visited = {}
        self.pos_traverser.extend_from_block(
            point.region_path_id, length, visited)
        visited = {node_id: min(self.graph.node_size(node_id), l)
                   for node_id, l in visited.items()}
        visited = {node_id: [0, l] for node_id, l in visited.items()}
        visited[point.region_path_id][0] = point.offset
        return visited

    def extend_interval_fast(self, interval, local_direction=1):
        """Direction: +1 forward, 0 both"""
        assert local_direction in [-1, 0, 1]
        interval.graph = self.graph
        extension_length = self.d - interval.length()
        areas = area_from_interval(interval, self.graph)
        end_position = interval.end_position
        new_areas = self.get_areas_from_point(
            end_position, extension_length)
        update_areas_fast(areas, new_areas)

        if local_direction == 0:
            end_position = interval.start_position
            new_areas = self.get_areas_from_point(
                end_position, self.d)
            update_areas_fast(areas, new_areas)
        return areas
