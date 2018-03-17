import logging
import pickle
from collections import deque

from offsetbasedgraph.graphtraverser import GraphTravserserBetweenNodes
import offsetbasedgraph as obg

from pyvg import Snarls


class SnarlGraph(obg.GraphWithReversals):

    def __init__(self, blocks, edges, id=None, parent=None,
                 children=[], start_node=None, end_node=None):
        super(SnarlGraph, self).__init__(blocks, edges)

        for child in children:
            assert child.id != id, "Child id %d = self.id %d" % (child.id, id)
        self._start_node = start_node
        self._end_node = end_node
        self.parent = parent
        self.children = children
        self.id = id
        self.create_children()

        self._length = None

    def get_next_nodes(self, node_id):
        if node_id not in self.adj_list:
            return []
        return self.adj_list[node_id]

    def get_previous_nodes(self, node_id):
        return self.reverse_adj_list[node_id]

    def node_size(self, node_id):
        length = self.blocks[abs(node_id)].length()
        assert length > 0
        return length

    def get_all_leaf_node_ids(self):
        leaf_nodes = [node_id for node_id, obj in self.blocks.items()
                      if isinstance(obj, obg.Block)]
        subgraphs = [obj for obj in self.blocks.values()
                     if isinstance(obj, SnarlGraph)]

        for subgraph in subgraphs:
            leaf_nodes.extend(subgraph.get_all_leaf_node_ids())
        return leaf_nodes

    def create_children(self):

        n_childs = len(self.children)
        for i, child in enumerate(self.children):
            if self.id == "top_level" and (i % 20000 == 0 or i < 4):
                logging.info("Snarl %s of %s" % (i, n_childs))
            assert child.id != self.id, "Child ID %d equal as parent" % child.id
            child_blocks, child_graph = SnarlGraph.create_from_simple_snarl(
                child, self)

            if len(child_blocks) == 0:
                continue

            self.blocks[child.id] = child_graph

            for node_id in child_blocks.keys():
                if node_id == child.start or node_id == child.end:
                    continue

                del self.blocks[abs(node_id)]
                self._delete_edges_to_and_from_node(node_id)
                self._delete_edge(child.start, node_id)
                self._delete_edge(child.start, -node_id)
                self._delete_edge(-child.end, node_id)
                self._delete_edge(-child.end, -node_id)

            # add new edges
            self.adj_list[child.start].append(child.id)
            self.reverse_adj_list[-child.id] = [-child.start]

            self.adj_list[child.id] = [child.end]
            self.reverse_adj_list[-child.end].append(-child.id)

    def _remove_edges(self, node_id):
        if node_id in self.adj_list:
            edges = self.adj_list[node_id]
            for edge in edges:
                self.reverse_adj_list[-edge].remove(-node_id)
            del self.adj_list[node_id]

        if -node_id in self.reverse_adj_list:
            edges = self.reverse_adj_list[-node_id]
            for edge in edges:
                self.adj_list[-edge].remove(node_id)
            del self.reverse_adj_list[-node_id]

    def _delete_edges_to_and_from_node(self, node_id):
        self._remove_edges(node_id)
        self._remove_edges(-node_id)
        return

    def _delete_edge(self, from_node, to_node):
        if to_node in self.adj_list[from_node]:
            self.adj_list[from_node].remove(to_node)

        if -from_node in self.reverse_adj_list[-to_node]:
            self.reverse_adj_list[-to_node].remove(-from_node)

    def _sanitize_graph(self):
        # Check that all edges goes to blocks that exists
        for block, adj_list in self.adj_list.items():
            for edge in adj_list:
                assert edge in self.blocks, "Block from %d going to %d. %d does not exist." % (block, edge, edge)

    @classmethod
    def create_from_simple_snarl(cls, simple_snarl, parent_snarl_graph):
        traverser = GraphTravserserBetweenNodes(parent_snarl_graph)
        subgraph = traverser.get_snarl_subgraph(
            simple_snarl.start, simple_snarl.end, include_start_and_end=True)

        # Remove start and end blocks (keep edges)
        del subgraph.blocks[simple_snarl.start]
        del subgraph.blocks[simple_snarl.end]

        return subgraph.blocks, SnarlGraph(
            subgraph.blocks, subgraph.adj_list, simple_snarl.id,
            parent=parent_snarl_graph,
            children=simple_snarl.children,
            start_node=simple_snarl.start,
            end_node=simple_snarl.end
        )

    def length(self):
        if self._length is not None:
            return self._length
        self._length = self._get_longest_path_length()
        return self._length

    def _get_longest_path_length(self):
        memo = self._create_path_length_dict()
        return memo[self._end_node]

    def _create_path_length_dict(self, forward=True):
        if forward:
            start_node = self._start_node
            end_node = self._end_node
            next_node_func = self.get_next_nodes
        else:
            start_node = -self._end_node
            end_node = -self._start_node
            next_node_func = self.get_previous_nodes

        stack = deque([(start_node, 0, 0)])
        memo = {}
        while stack:
            node_id, dist, n = stack.pop()
            next_nodes = next_node_func(node_id)
            for next_node in next_nodes:
                if next_node in memo and memo[next_node] >= dist:
                    continue
                memo[next_node] = dist
                if next_node == end_node:
                    continue
                new_dist = dist + self.node_size(next_node)
                stack.append((next_node, new_dist, n+1))

        return memo

    def get_linear_node_intervals(self):
        linear_node_intervals = {}
        forward_length_dict = self._create_path_length_dict()
        back_length_dict = self._create_path_length_dict(False)
        for node_id in self.blocks:
            if (node_id in forward_length_dict) and (-node_id in back_length_dict):
                start_length = forward_length_dict[node_id]
                end_length = back_length_dict[-node_id]
            elif (-node_id in forward_length_dict) and (node_id in back_length_dict):
                # Use reverse node
                start_length = forward_length_dict[-node_id]
                end_length = back_length_dict[node_id]
            else:
                logging.error(node_id)
                logging.error(forward_length_dict)
                logging.error(back_length_dict)
                raise Exception()

            path_length = start_length + end_length + self.node_size(node_id)
            scale_factor = self.length() / path_length
            if scale_factor < 1:
                raise Exception()

            linear_node_intervals[node_id] = (
                start_length * scale_factor,
                (start_length+self.node_size(node_id))*scale_factor)

        return linear_node_intervals

    def get_distance_dicts(self):
        starts_dict, ends_dict = {}, {}
        linear_node_intervals = self.get_linear_node_intervals()
        for node_id, block in self.blocks.items():
            # TODO: Also get reverse distances
            start, end = linear_node_intervals[node_id]
            if isinstance(block, obg.Block):
                starts_dict[node_id] = start
                ends_dict[node_id] = end
            else:
                sub_starts_dict, sub_ends_dict = block.get_distance_dicts()
                starts_dict.update(
                    {n_id: start+dist
                     for n_id, dist in sub_starts_dict.items()})
                ends_dict.update({n_id: start + dist
                                  for n_id, dist in sub_ends_dict.items()})
        for node_id in starts_dict.keys():
            assert starts_dict[node_id] != ends_dict[node_id]
        return starts_dict, ends_dict

    def find_node_ids(self, offset):
        sub_nodes = [node_id for node_id, interval
                     in self.linear_node_intervals.items()
                     if interval[0] < offset and interval[1] > offset]
        positions = []
        for node_id in sub_nodes:
            if isinstance(self.blocks[node_id],
                          obg.Block):
                positions.append(node_id)
            else:
                positions.extend(self.blocks[node_id].find_node_ids())

        return positions

    def to_file(self, file_name):
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def from_file(cls, file_name):
        with open("%s" % file_name, "rb") as f:
            obj = pickle.loads(f.read())
            assert isinstance(obj, cls)
            return obj


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

    def add_child(self, child):
        assert child.parent == self.id
        assert child.id != self.id
        self.children.append(child)

    def set_parent(self, parent_id):
        self.parent = parent_id

    def __str__(self):
        return "Snarl(%d, %d, id=%d, parent=%s)" % (
            self.start, self.end, self.id, self.parent)

    def __repr__(self):
        return self.__str__()

    def sanitize(self):
        for child in self.children:
            assert child.id != self.id, " Child %d has same id as self %d" % (child.id, self.id)


class SnarlGraphBuilder:

    def __init__(self, graph, snarls, id_counter=0):
        assert isinstance(snarls, dict), \
            "Should be dict with ids as keys and SimpleSnarls as values"
        assert id_counter > 0,\
            "ID counter must be set to higher than highest ID in graph"
        self.graph = graph
        self.snarls = snarls  # dict of snarls, id: SimpleSnarls
        self.id_counter = id_counter

    def build_snarl_graphs(self):
        logging.info("Finding top level snarls")
        top_level_snarls = [snarl for snarl in self.snarls.values()
                            if snarl.parent is None]

        # Find all snarls with 0 nodes
        logging.info("Finding snarls with 0 nodes")
        n_zero_nodes = 0
        for snarl in self.snarls.values():
            start = snarl.start
            end = snarl.end

            if end in self.graph.adj_list[start]:
                n_zero_nodes += 1

        # Add dummy start and end to graph
        new_start = self.id_counter
        new_end = new_start + 1

        for block in self.graph.get_first_blocks():
            self.graph._add_edge(new_start, block)

        for block in self.graph.get_last_blocks():
            self.graph._add_edge(block, new_end)

        logging.info("%d snarls with 0 nodes" % n_zero_nodes)
        logging.info("%d top level" % len(top_level_snarls))
        logging.info("%d snarls total" % len(self.snarls))

        whole_graph_snarl = SnarlGraph(self.graph.blocks,
                                       self.graph.adj_list,
                                       "top_level", parent=None,
                                       children=top_level_snarls,
                                       start_node=new_start,
                                       end_node=new_end)
        return whole_graph_snarl

    @classmethod
    def from_vg_snarls(cls, graph, vg_snarls_file_name):
        logging.info("Building snarlgraph from vg graph")
        snarls = Snarls.from_vg_snarls_file(vg_snarls_file_name).snarls
        start_end_mapping = {}
        simple_snarls = {}

        max_graph_id = max([id for id in graph.blocks.keys()])
        id_counter = max_graph_id + 1
        i = 0
        for snarl in snarls:

            assert snarl.start.backward == snarl.end.backward
            start = snarl.start.node_id
            end = snarl.end.node_id

            if snarl.start.backward:
                start = snarl.end.node_id
                end = snarl.start.node_id

            simple_snarls[id_counter] = SimpleSnarl(
                start, end, id=id_counter, parent=snarl.parent)

            start_end_mapping["%d-%d" % (start, end)] = id_counter
            id_counter += 1

            if i % 500000 == 0:
                logging.info("Reading snarl %d from vg" % i)
            i += 1

        logging.info("Setting parent snarls")
        # Set parent
        for id, snarl in simple_snarls.items():
            if snarl.parent is not None:
                parent = snarl.parent
                if parent.start.node_id == 0:
                    snarl.set_parent(None)
                else:
                    mapping_id = "%d-%d" % (
                        parent.start.node_id, parent.end.node_id)
                    if mapping_id in start_end_mapping:
                        parent_id = start_end_mapping[mapping_id]
                        parent = simple_snarls[parent_id]
                        snarl.set_parent(parent.id)
                        assert snarl.parent == parent.id
                        parent.add_child(snarl)
                        assert snarl in parent.children
                        assert parent.id != snarl.id

                    else:
                        logging.warning("Found with parent that is missing")
                        logging.warning(snarl)

        logging.info("Done reading snarls")
        return cls(graph, simple_snarls, id_counter)
