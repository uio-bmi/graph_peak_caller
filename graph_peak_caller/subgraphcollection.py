from .sparsepileup import SparsePileup
from .extender import Areas
import numpy as np

class ConnectedAreas(Areas):
    def __init__(self, graph, areas=None):
        super(ConnectedAreas, self).__init__(graph, areas)

    def touches_node(self, other_node):
        if other_node in self.areas:
            return True
            #raise Exception("Checking whether connected areas touching a node id it already contains")

        for node, starts_and_ends in self.areas.items():
            if starts_and_ends[-1] == self.graph.node_size(node):
                if other_node in self.graph.adj_list[node]:
                    return True
                elif -node in self.graph.reverse_adj_list[-other_node]:
                    return True
            if starts_and_ends[0] == 0:
                if node in self.graph.adj_list[other_node]:
                    return True
                elif -other_node in self.graph.reverse_adj_list[-node]:
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

        collection = cls.__init__(graph)

        areas = pileup.find_valued_areas()
        for node_id, starts_and_ends in areas.items():
            for i in range(0, len(starts_and_ends) // 2):
                collection.add_area(node_id,
                                    starts_and_ends[i*2],
                                    starts_and_ends[i*2+1])

    def _subgraphs_touching_node(self, node_id):
        touching_subgraphs = []
        for subgraph in self.subgraphs:
            if subgraph.touches_node(node_id):
                touching_subgraphs.append(subgraph)

        return touching_subgraphs

    def add_area(self, node_id, start, end):
        assert node_id in self.graph.blocks
        assert start >= 0 and start < end
        assert end <= self.graph.node_size(node_id)

        touching_subgraphs = self._subgraphs_touching_node(node_id)
        #print("Adding %d, %d, %d" % (node_id, start, end))
        #print("Touching subgraphs")
        #print(touching_subgraphs)

        if len(touching_subgraphs) == 0:
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






