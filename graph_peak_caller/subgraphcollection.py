from .extender import Areas
import numpy as np
from offsetbasedgraph import Position

class ConnectedAreas(Areas):
    def __init__(self, graph, areas=None):
        super(ConnectedAreas, self).__init__(graph, areas)

    def touches_area(self, other_node, start, end):
        graph = self.graph

        if start > 0 and end < graph.node_size(other_node):
            if other_node not in self.areas:
                return False

        if start == 0:
            if not self._is_start_position(other_node, start):
               return True

        if end == self.graph.node_size(other_node):
            if not self._is_end_position(other_node, end):
                return True

        return False

    def __add__(self, other):
        # Adds another connected area to this one
        for node_id, starts_and_ends in other.areas.items():
            self.add_areas_for_node(node_id, starts_and_ends)
        return self

class SubgraphCollection(object):

    def __init__(self, graph, subgraphs=None):
        self.graph = graph
        # Every subgraph is a list of connected areas
        if subgraphs is not None:
            self.subgraphs = subgraphs
        else:
            self.subgraphs = []

    def remove(self, subgraph):
        self.subgraphs.remove(subgraph)

    @classmethod
    def from_pileup(cls, graph, pileup):
        collection = cls(graph)
        areas = pileup.find_valued_areas(1)
        for node_id, starts_and_ends in areas.items():
            for i in range(0, len(starts_and_ends) // 2):
                collection.add_area(node_id,
                                    starts_and_ends[i*2],
                                    starts_and_ends[i*2+1])
        return collection

    def _subgraphs_touching_area(self, node_id, start, end):
        touching_subgraphs = []
        for subgraph in self.subgraphs:
            if subgraph.touches_area(node_id, start, end):
                touching_subgraphs.append(subgraph)

        return touching_subgraphs

    def add_area(self, node_id, start, end):
        assert node_id in self.graph.blocks
        assert start >= 0 and start < end
        assert end <= self.graph.node_size(node_id)

        touching_subgraphs = self._subgraphs_touching_area(node_id, start, end)

        if len(touching_subgraphs) == 0:
            print("Creating new subgraph for %d, %d, %d" % (node_id, start, end))
            new_subgraph = ConnectedAreas(self.graph, {node_id: np.array([start, end])})
            self.subgraphs.append(new_subgraph)
        elif len(touching_subgraphs) == 1:
            touching_subgraphs[0].add_areas_for_node(node_id, np.array([start, end]))
        elif len(touching_subgraphs) > 1:
            # Merge all touching subgraphs, then add the area
            new_subgraph = touching_subgraphs[0]
            for touching_subgraph in touching_subgraphs[1:]:
                new_subgraph = new_subgraph + touching_subgraph
                self.remove(touching_subgraph)

            new_subgraph.add_areas_for_node(node_id, np.array([start, end]))

    def to_file(self, file_name):
        f = open(file_name, "w")
        lines = (subgraph.to_file_line() for subgraph in self.subgraphs)
        f.writelines(lines)
        f.close()
