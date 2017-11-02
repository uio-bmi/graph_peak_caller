from offsetbasedgraph import Graph
from pyvg import Snarls
import logging
from filecache import filecache
from offsetbasedgraph.graphtraverser import GraphTravserserBetweenNodes


class SnarlGraph(Graph):

    def __init__(self, blocks, edges, id=None, parent=None, children=[]):
        super(SnarlGraph, self).__init__(blocks, edges)

        for child in children:
            assert child.id != id, "Child id %d = self.id %d" % (child.id, id)

        self.parent = parent
        self.children = children  # Simple snarls (not yet created as SnarlGraphs)
        #print("%d children" % len(self.children))
        self.id = id
        self.create_children()

    def create_children(self):
        #print("Creating for children from %s. %d children" % (self.id, len(self.children)))

        for child in self.children:
            #print("    Child: %d, %s" % (child.id, child))

            assert child.id != self.id, "Child ID %d equal as parent" % child.id

            child_blocks, child_graph = SnarlGraph.create_from_simple_snarl(child, self)
            #print("         %d blocks in child" % len(child_blocks))
            if len(child_blocks) == 0:
                #print("      0 child blocks. Continue")
                continue

            self.blocks[child.id] = child_graph
            #print("  added block %d" % child.id)

            #print("             Child subgraph from %s:" % child)
            #print(child_graph)

            for node_id in child_blocks.keys():

                #assert node_id != 598821

                #print("         Removing block %d" % node_id)
                if node_id == child.start or node_id == child.end:
                    continue

                del self.blocks[abs(node_id)]

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

                if -node_id in self.adj_list[child.start]:
                    self.adj_list[child.start].remove(-node_id)

                if -node_id in self.reverse_adj_list[-child.end]:
                    self.reverse_adj_list[-child.end].remove(-node_id)

                if node_id in self.reverse_adj_list[-child.end]:
                    self.reverse_adj_list[-child.end].remove(node_id)

            # add new edges
            self.adj_list[child.start].append(child.id)
            self.reverse_adj_list[-child.id] = [-child.start]

            #print("     adding edge  %d to %d" %  (child.start, child.id))
            #print("     adding reversee edge  %d to %d" %  (-child.id, -child.start))

            self.adj_list[child.id] = [child.end]
            self.reverse_adj_list[-child.end].append(-child.id)


            #print("     adding edge  %d to %d" %  (child.id, child.end))
            #print("     adding reversee edge  %d to %d" %  (-child.end, -child.id))

            #self._sanitize_graph()

    def _sanitize_graph(self):
        # Check that all edges goes to blocks that exists
        for block, adj_list in self.adj_list.items():
            for edge in adj_list:
                assert edge in self.blocks, "Block from %d going to %d. %d does not exist." % (block, edge, edge)

    @classmethod
    def create_from_simple_snarl(cls, simple_snarl, parent_snarl_graph):
        traverser = GraphTravserserBetweenNodes(parent_snarl_graph)
        print_debug = False
        #subgraph = traverser.get_greedy_subgraph_between_nodes(simple_snarl.start, simple_snarl.end, include_start_and_end=False, print_debug=print_debug)
        subgraph = traverser.get_snarl_subgraph(simple_snarl.start, simple_snarl.end, include_start_and_end=False)

        simple_snarl.sanitize()

        return subgraph.blocks, SnarlGraph(subgraph.blocks, subgraph.adj_list, simple_snarl.id, parent=parent_snarl_graph, children=simple_snarl.children)

    def length(self):
        pass



class SimpleSnarl():
    def __init__(self, start, end, id, parent=None, children=None):
        self.start = start
        self.end = end
        self.id = id
        self.parent = parent
        if children is None:
            self.children = []
        else:
            self.children = children

        assert self.parent != self.id

        #self.sanitize()

    def add_child(self, child):
        #print("    Childs children 1: %s" % child.children)
        assert child.parent == self.id
        assert child.id != self.id
        #print("    Childs children 2: %s" % child.children)
        self.children.append(child)
        #print("    Childs children 3: %s" % child.children)
        #print("   Adding child %d to parent  %d" % (child.id, self.id))
        #print("    %s" % self.children)
        #print("    Childs children: %s" % child.children)
        child.sanitize()

    def set_parent(self, parent_id):
        self.parent = parent_id

    def __str__(self):
        return "Snarl(%d, %d, id=%d, parent=%s)" % (self.start, self.end, self.id, self.parent)

    def __repr__(self):
        return self.__str__()

    def sanitize(self):
        #print("Sanitizing %d." % self.id)
        for child in self.children:
            #print("   Checking child %d" % child.id)
            assert child.id != self.id, " Child %d has same id as self %d" % (child.id, self.id)



class SnarlGraphBuilder:

    def __init__(self, graph, snarls):
        self.graph = graph
        self.snarls = snarls  # dict of snarls, id: SimpleSnarls

    def build_snarl_graphs(self):
        top_level_snarls = [snarl for snarl in self.snarls.values() if snarl.parent is None]

        # Find all snarls with 0 nodes
        n_zero_nodes = 0
        for snarl in self.snarls.values():
            start = snarl.start
            end = snarl.end

            if end in self.graph.adj_list[start]:
                n_zero_nodes += 1

        print("%d snarls with 0 nodes" % n_zero_nodes)

        print("%d top level" % len(top_level_snarls))
        print("%d snarls total" % len(self.snarls))
        #return
        whole_graph_snarl = SnarlGraph(self.graph.blocks, self.graph.adj_list, "top_level", None, top_level_snarls)
        return whole_graph_snarl

    @classmethod
    #@filecache(24*60*60)
    def from_vg_snarls(cls, graph, vg_snarls_file_name):

        snarls = Snarls.from_vg_snarls_file(vg_snarls_file_name).snarls
        start_end_mapping = {}

        #for snarl in snarls:
        #    assert len(graph.adj_list[snarl.start.node_id]) == 1, \
        #        "Start node %d has multiple edges out: %s" % (snarl.start.node_id, graph.adj_list[snarl.start.node_id])

        simple_snarls = {}
        max_graph_id = max([id for id in graph.blocks.keys()])
        id_counter = max_graph_id + 1
        for snarl in snarls:

            assert snarl.start.backward == snarl.end.backward

            start = snarl.start.node_id
            end = snarl.end.node_id

            if snarl.start.backward:
                #print("Reversed snarl %d, %d" % (start, end))
                start = snarl.end.node_id
                end = snarl.start.node_id

            simple_snarls[id_counter] = SimpleSnarl(start, end, id=id_counter, parent=snarl.parent)

            start_end_mapping["%d-%d" % (start, end)] = id_counter
            id_counter += 1

        # Set parent
        for id, snarl in simple_snarls.items():

            #if snarl.start == 485840:
            #    print(snarl)
            #    return

            if snarl.parent is not None:
                parent = snarl.parent
                if parent.start.node_id == 0:
                    #print("No parent")
                    snarl.set_parent(None)
                else:
                    mapping_id = "%d-%d" % (parent.start.node_id, parent.end.node_id)
                    if mapping_id in start_end_mapping:
                        parent_id = start_end_mapping[mapping_id]
                        parent = simple_snarls[parent_id]

                        #print("Creating child from %d to %d" % (parent.id, snarl.id))
                        #print("N children in snarl: %d" % len(snarl.children))
                        #snarl.sanitize()
                        snarl.set_parent(parent.id)
                        #snarl.sanitize()
                        assert snarl.parent == parent.id
                        #print("Snarls children: %s " % snarl.children)
                        #print("Parent: %s" % parent)
                        parent.add_child(snarl)
                        assert snarl in parent.children

                        snarl.sanitize()
                        #simple_snarls[parent].sanitize()
                        assert parent.id != snarl.id

                    else:
                        logging.warning("Found with parent that is missing")
                        logging.warning(snarl)

        #for snarl in simple_snarls.values():
        #    snarl.sanitize()

        return cls(graph, simple_snarls)
