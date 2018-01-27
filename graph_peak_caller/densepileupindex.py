from offsetbasedgraph.graphtraverser import GraphTraverser
from collections import defaultdict

class GraphIndex(object):
    def __init__(self, index):
        self.index = index

    def create_from_graph(self, graph, length=100):
        self._create_index(graph, length)

    def get_node(self, node):
        yield (node, 0)
        for child in self.index[node]:
            yield child

    def _create_index(self, graph, length):
        self.index = defaultdict(list)
        for node in graph:
            self._traverse_from_node(graph, node, length)
            self._traverse_from_node(graph, -node, length)

    def _traverse_from_node(self, graph, node, length):
        if node > 0:
            traverser = GraphTraverser(graph, direction=1)
        else:
            traverser = GraphTraverser(graph, direction=-1)

        visited = defaultdict(int)
        for next_node in traverser.adj_list[node]:
            traverser.extend_from_block(next_node, length, visited)

        for to_node, to_length in visited:
            self.index[node].append((to_node, to_length))

    @classmethod
    def from_file(cls, file_name):
        pass

    def to_file(self, file_name):
        pass


class GraphExtender(object):
    def __init__(self, graph_index):
        self.index = graph_index

    def extend_from_position(self, node, offset, extension_length):
        return ((node, length-offset) for node, length in
                self.index.get_node(node) if length-offset < extension_length)


class DensePileupExtender(object):
    def __init__(self, graph_extender, dense_pileup):
        self.graph_extender = graph_extender
        self.pileup = dense_pileup

    def extend_from_position(self, node, offset, extension_length):
        extensions = self.graph_extender.extend_from_position(node, offset, extension_length)
        for extension in extensions:
            print("Extension: %d,%d" % (extension))
            to_node, length_to_start = extension
            to_node_end = extension_length - length_to_start
            to_node_start = max(0, -length_to_start)
            print("start, end: %d, %d" % (to_node_start, to_node_end))
            yield self.pileup.data.node_range_to_value_indexes(to_node, to_node_start, to_node_end)



"""

TODO:
    - Create graph index from graph
    - Create sample extender

For every node:

    List of nodes with start less than length away from node start
        For every node in this list:
            Contain start and end position in dense pileup values


    Given position in graph, will ask for all possible positions in dense pileup that are n basepairs (extension_size) away
    Position: node, offset
        1) Ask node for list of children
        2) For every node in children
            1) Get start, end densepileup indexes
            2) Adjust start end to:
                n_start = max(end, start + offset)
                n_end = min(n_start + extension_size, end)
"""



