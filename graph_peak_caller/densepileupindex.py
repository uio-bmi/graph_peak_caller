from offsetbasedgraph.graphtraverser import LengthGraphTraverser
from collections import defaultdict
import pickle
import logging

class GraphIndex(object):
    def __init__(self, index):
        self.index = index

    @classmethod
    def create_from_graph(cls, graph, length=100):
        index = GraphIndex(None)
        max_node_size = graph.max_node_size()
        assert length > max_node_size, \
            "Length should be much larger than max node size (%d)" % max_node_size
        index._create_index(graph, length)
        return index

    def get_node(self, node):
        yield (node, 0)
        for child in self.index[node]:
            yield child

    def _create_index(self, graph, length):
        self.index = defaultdict(list)
        i = 0
        for node in graph.blocks:
            if i % 20000 == 0:
                logging.info("Processing node %d of %d" % (i, len(list(graph.blocks.keys()))))
            i += 1
            self._traverse_from_node(graph, node, length)
            self._traverse_from_node(graph, -node, length)

    def _traverse_from_node(self, graph, node, length):
        #print("Traversing from %d" % node)
        if node > 0:
            traverser = LengthGraphTraverser(graph, length, direction=1)
        else:
            traverser = LengthGraphTraverser(graph, length, direction=-1)

        distances = traverser.extend_from_block(node)

        for to_node, to_length in distances.items():
            self.index[node].append((to_node, to_length))

    @classmethod
    def from_file(cls, file_base_name):
        with open(file_base_name + ".obgindex", "rb") as f:
            index = pickle.load(f)
        return GraphIndex(index)

    def to_file(self, file_base_name):
        with open(file_base_name + ".obgindex", "wb") as f:
            pickle.dump(self.index, f)

    def __str__(self):
        out = "GraphIndex:\n"
        for node in self.index.keys():
            out += "   %d: %s\n" %\
                   (node, ', '.join([str(dist) for dist in self.get_node(node)]))
        return out

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        for node in self.index.keys():
            if list(self.get_node(node)) != list(other.get_node(node)):
                return False
        for node in other.index.keys():
            if list(self.get_node(node)) != list(other.get_node(node)):
                return False
        return True


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

    def extend_from_position(self, node, offset, extension_length, touched_nodes=None):
        extensions = self.graph_extender.extend_from_position(node, offset, extension_length)
        for extension in extensions:
            #print("Extension: %d,%d" % (extension))
            to_node, length_to_start = extension
            to_node_end = extension_length - length_to_start
            to_node_start = max(0, -length_to_start)
            #print("start, end: %d, %d" % (to_node_start, to_node_end))
            if touched_nodes is not None:
                touched_nodes.add(abs(to_node))
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



