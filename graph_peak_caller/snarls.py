from offsetbasedgraph import Graph
from pyvg import Snarls
import logging
from filecache import filecache
from offsetbasedgraph.graphtraverser import GraphTravserserBetweenNodes


class SnarlGraph(Graph):

    def __init__(self, blocks, edges, id=None, parent=None, children=[]):
        super(SnarlGraph, self).__init__(blocks, edges)

        self.parent = parent
        self.children = children  # Simple snarls (not yet created as SnarlGraphs)
        #print("%d children" % len(self.children))
        self.id = id
        self.create_children()

    def create_children(self):
        print("Creating for children from %s" % self.id)
        for child in self.children:
            print("    Child: %d, %s" % (child.id, child))

            child_graph = SnarlGraph.create_from_simple_snarl(child, self)
            print("      created child_graph")
            self.blocks[child.id] = child_graph

            print("        Found subgraph:")
            print(child_graph)

            for node_id in child_graph.blocks.keys():

                if node_id == child.start or node_id == child.end:
                    continue

                del self.blocks[node_id]

                if node_id in self.adj_list:
                    del self.adj_list[node_id]

                if -node_id in self.adj_list:
                    del self.adj_list[-node_id]

                if node_id in self.reverse_adj_list:
                    del self.reverse_adj_list[node_id]

                if -node_id in self.reverse_adj_list:
                    del self.reverse_adj_list[-node_id]

                if node_id in self.adj_list[child.start]:
                    self.adj_list[child.start].remove(node_id)

                if -node_id in self.reverse_adj_list[-child.end]:
                    self.reverse_adj_list[-child.end].remove(-node_id)

            # add new edges
            self.adj_list[child.start].append(child.id)
            self.reverse_adj_list[-child.id] = [-child.start]

            self.adj_list[child.id] = [child.end]
            self.reverse_adj_list[-child.end].append(-child.id)

    @classmethod
    def create_from_simple_snarl(cls, simple_snarl, parent_snarl_graph):
        traverser = GraphTravserserBetweenNodes(parent_snarl_graph)
        subgraph = traverser.get_greedy_subgraph_between_nodes(simple_snarl.start, simple_snarl.end, include_start_and_end=False)
        return SnarlGraph(subgraph.blocks, subgraph.adj_list, simple_snarl.id, parent_snarl_graph, simple_snarl.children)

    def length(self):
        pass


class SimpleSnarl():
    def __init__(self, start, end, id, parent=None, children=[]):
        self.start = start
        self.end = end
        self.id = id
        self.parent = parent
        self.children = children

    def add_child(self, child):
        self.children.append(child)

    def set_parent(self, parent_id):
        self.parent = parent_id

    def __str__(self):
        return "Snarl(%d, %d, id=%d, parent=%s)" % (self.start, self.end, self.id, self.parent)

    def __repr__(self):
        return self.__str__()


class SnarlGraphBuilder:

    def __init__(self, graph, snarls):
        self.graph = graph
        self.snarls = snarls  # dict of snarls, id: SimpleSnarls

    def build_snarl_graphs(self):
        top_level_snarls = [snarl for snarl in self.snarls.values() if snarl.parent is None]
        print("%d top level" % len(top_level_snarls))
        whole_graph_snarl = SnarlGraph(self.graph.blocks, self.graph.adj_list, "top_level", None, top_level_snarls)
        return whole_graph_snarl

    @classmethod
    #@filecache(24*60*60)
    def from_vg_snarls(cls, graph, vg_snarls_file_name):

        snarls = Snarls.from_vg_snarls_file(vg_snarls_file_name).snarls
        start_end_mapping = {}

        simple_snarls = {}
        max_graph_id = max([id for id in graph.blocks.keys()])
        id_counter = max_graph_id
        for snarl in snarls:
            simple_snarls[id_counter] = SimpleSnarl(snarl.start.node_id, snarl.end.node_id, id_counter)
            start_end_mapping["%d-%d" % (snarl.start.node_id, snarl.end.node_id)] = id_counter
            id_counter += 1

        # Set parent
        for id, snarl in simple_snarls.items():
            if snarl.parent is not None:
                parent = snarl.parent
                if parent.start.node_id == 0:
                    pass
                else:
                    mapping_id = "%d-%d" % (parent.start.node_id, parent.end.node_id)
                    if mapping_id in start_end_mapping:
                        parent_id = start_end_mapping[mapping_id]
                        simple_snarls[id].set_parent(parent_id)
                        simple_snarls[parent_id].add_children(simple_snarls[id])
                        #print("Parent: %s" % mapping_id)
                    else:
                        logging.warning("Found with parent that is missing")
                        logging.warning(snarl)

        return cls(graph, simple_snarls)
