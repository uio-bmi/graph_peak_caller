from .extender import Areas
import numpy as np
import pickle
import logging
from collections import OrderedDict
from .areas import BinaryContinousAreas
from .extender import Areas
import itertools
from collections import defaultdict

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

    def contains_interval(self, interval):
        intervals = self.to_simple_intervals()
        overlap = 0
        for internal_interval in intervals:
            overlap += interval.overlap(internal_interval)

        assert overlap <= interval.length()

        if overlap == interval.length():
            return True

        return False

    def n_basepairs(self):
        intervals = self.to_simple_intervals()
        n = 0
        for interval in intervals:
            n += interval.length()
        return n


class BCConnectedAreas(BinaryContinousAreas):
    def __init__(self, graph, id, areas=None):
        super(BCConnectedAreas, self).__init__(graph)
        self.id = id
        for node, startend in areas.items():
            self.add(node, startend[0], startend[1])

    def __hash__(self):
        return self.id

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
        self.merge_with_other(other)
        return self

    def contains_interval(self, interval):
        intervals = self.to_simple_intervals()
        overlap = 0
        for internal_interval in intervals:
            overlap += interval.overlap(internal_interval)

        assert overlap <= interval.length()

        if overlap == interval.length():
            return True

        return False

    def n_basepairs(self):
        intervals = self.to_simple_intervals()
        n = 0
        for interval in intervals:
            n += interval.length()
        return n


class SubgraphCollection(object):

    def __init__(self, graph, subgraphs=None):
        self.graph = graph
        self._id_counter = 0
        if subgraphs is not None:
            self.subgraphs = subgraphs
        else:
            self.subgraphs = []

        self._node_end_index = defaultdict(set)  # Index from node id to subgraphs touching end
        self._node_start_index = defaultdict(set)

    def remove(self, subgraph):
        self.subgraphs.remove(subgraph)

    def add(self, subgraph):
        pass

    @classmethod
    def from_pileup(cls, graph, pileup):
        collection = cls(graph)
        logging.info("Finding valued areas")
        areas = pileup.find_valued_areas(1)
        n = 0
        for node_id, starts_and_ends in areas.items():
            if n % 50000 == 0:
                logging.info("Adding node %d to subgraph" % n)
            n += 1
            for i in range(0, len(starts_and_ends) // 2):
                collection.add_area(node_id,
                                    starts_and_ends[i*2],
                                    starts_and_ends[i*2+1])
        return collection

    def _subgraphs_touching_area(self, node_id, start, end):
        out = set()
        if start == 0:
            #print("   Checking start of node %d" % node_id)
            out = self._node_start_index[node_id]
            #print("   Found:")
            #print(self._node_start_index[node_id])
        if end == self.graph.node_size(node_id):
            out.update(self._node_end_index[node_id])
            #print("  Checking end of node %d" % node_id)
            #print("   Found:")
            #print(self._node_end_index[node_id])
        #print("   Found in total %d touching subgraphs" % len(out))
        #return list(set(out))
        return out

        touching_subgraphs = []
        for subgraph in self.subgraphs:
            if subgraph.touches_area(node_id, start, end):
                touching_subgraphs.append(subgraph)
        return touching_subgraphs

    def __iter__(self):
        return iter(self.subgraphs)

    def add_indexes(self, node, start, end, subgraph):
        # Checks if added node id,start,end touches nodes. Link index to subgraph
        if start == 0:
            for in_node in self.graph.adj_list[-node]:
                #print("    Adding subgraph to end index of %d and start of %d" % (in_node, -in_node))
                self._node_end_index[in_node].add(subgraph)
                self._node_start_index[-in_node].add(subgraph)
            for in_node in self.graph.reverse_adj_list[-node]:
                #print("    Adding subgraph to end index of %d and end of %d" % (-in_node, in_node))
                self._node_end_index[-in_node].add(subgraph)
                self._node_start_index[in_node].add(subgraph)

        if end == self.graph.node_size(node):
            for in_node in self.graph.adj_list[node]:
                #print("    Adding subgraph to start index of %d and end index of %d" % (in_node, -in_node))
                self._node_start_index[in_node].add(subgraph)
                self._node_end_index[-in_node].add(subgraph)
            for in_node in self.graph.reverse_adj_list[node]:
                #print("    Adding subgraph to end index of %d and start of %d" % (-in_node, in_node))
                self._node_end_index[-in_node].add(subgraph)
                self._node_start_index[in_node].add(subgraph)



    def add_area(self, node_id, start, end):
        assert node_id in self.graph.blocks
        assert start >= 0 and start < end
        assert end <= self.graph.node_size(node_id)
        #print("Adding node %d, startend %d %d" % (node_id, start, end))
        #print("%d subgraphs so far" % len(self.subgraphs))
        if start > 0 and end < self.graph.node_size(node_id):
            # Internal, will never touch anything
            touching_subgraphs = []
        else:
            touching_subgraphs = list(self._subgraphs_touching_area(node_id, start, end))

        self._id_counter += 1
        if len(touching_subgraphs) == 0:
            new_subgraph = BCConnectedAreas(
                self.graph,
                self._id_counter,
                {node_id: np.array([start, end])})
            self.subgraphs.append(new_subgraph)
            self.add_indexes(node_id, start, end, new_subgraph)

        elif len(touching_subgraphs) == 1:
            touching_subgraphs[0].add(
                node_id, start, end)
            #print("Subgraph after adding")
            #print(touching_subgraphs[0])
            self.add_indexes(node_id, start, end, touching_subgraphs[0])

        elif len(touching_subgraphs) > 1:
            #print("  Found multiple subgraphs")
            # Merge all touching subgraphs, then add the area
            new_subgraph = touching_subgraphs[0]
            for touching_subgraph in touching_subgraphs[1:]:
                new_subgraph = new_subgraph + touching_subgraph
                if touching_subgraph in self.subgraphs:
                    self.remove(touching_subgraph)

            new_subgraph.add(node_id, start, end)
            self.add_indexes(node_id, start, end, new_subgraph)

    def to_file(self, file_name):
        f = open(file_name, "w")
        lines = (subgraph.to_file_line() for subgraph in self.subgraphs)
        f.writelines(lines)
        f.close()

    def contains_interval(self, interval):
        for connected_area in self.subgraphs:
            if connected_area.contains_interval(interval):
                return True

        return False

    @classmethod
    def from_pickle(cls, file_name, graph=None):
        with open("%s" % file_name, "rb") as f:
            obj = pickle.loads(f.read())
            assert isinstance(obj, cls)

            for subgraph in obj.subgraphs:
                subgraph.graph = graph

            return obj

    def to_pickle(self, file_name):
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)


class SubgraphCollectionPartiallyOrderedGraph(SubgraphCollection):
    @classmethod
    def create_from_pileup(cls, graph, pileup):
        builder = SubgraphCollectionPartiallyOrderedGraphBuilder(
                    graph, pileup)
        return builder.build()


class SingleArea:
    def __init__(self, node, start, end, is_start, is_end):
        self.node = node
        self.start = start
        self.end = end
        self.is_start = is_start
        self.is_end = is_end

    def __str__(self):
        return "Node %d, start/end: %d/%d. Is start: %s, is end:%s" % \
               (self.node, self.start, self.end, self.is_start, self.is_end)

    def __repr__(self):
        return self.__str__()


class SubgraphCollectionPartiallyOrderedGraphBuilder():

    def __init__(self, graph, pileup):
        self.graph = graph
        self.pileup = pileup
        from .densepileup import DensePileupData
        assert isinstance(self.pileup.data, DensePileupData)
        areas = pileup.find_valued_areas(1)
        self.areas = Areas(self.graph, areas)
        self.single_areas = []

    def _make_single_areas(self):
        logging.info("Making single areas")
        sorted_nodes = sorted(self.areas.areas.keys())
        n = 0
        for node_id in sorted_nodes: #, starts_and_ends in self.areas.areas.items():
            #print("Checking node %d" % node_id)
            starts_and_ends = self.areas.areas[node_id]
            if n % 1000000 == 0:
                logging.info("Adding node %d to subgraph" % n)
            n += 1
            for i in range(0, len(starts_and_ends) // 2):
                start = starts_and_ends[i*2]
                end = starts_and_ends[i*2+1]
                is_start = self.areas._is_start_position(node_id, start)
                #print("Checking end %d, %d" % (node_id, end))
                is_end = self.areas._is_end_position(node_id, end)
                yield SingleArea(node_id, start, end, is_start, is_end)

    def _get_connected_subgraphs(self, subgraphs, area):
        connected = []
        #print("Getting connected subgraphs")
        for subgraph in subgraphs:
            #print("  Checking subgraph")
            if not area.is_start:
                for prev_node in self.graph.reverse_adj_list[-area.node]:
                    #print("        Checking prev node %d" % prev_node)
                    if -prev_node in subgraph.full_areas:
                        connected.append(subgraph)
                    if prev_node in subgraph.starts:
                        connected.append(subgraph)
            elif not area.is_end:
                for next_node in self.graph.adj_list[area.node]:
                    #print("        Checking next node %d" % next_node)
                    if next_node in subgraph.starts:
                        connected.append(subgraph)
                    if next_node in subgraph.full_areas:
                        connected.append(subgraph)

        return connected

    def _create_subgraphs(self):
        from collections import deque
        active_subgraphs = deque()

        logging.info("Creating subgraphs after single areas are created")
        peaks = []
        current_peak = BinaryContinousAreas(self.graph)
        prev_was_end = False
        prev_node = 0
        for single_area in self.single_areas:
            #print("Single area: %s " % single_area)

            connected = self._get_connected_subgraphs(active_subgraphs, single_area)
            if len(connected) > 1:
                merged = connected[0]
                for subgraph in connected[1:]:
                    merged.merge_with_other(subgraph)

                connected = [merged]

            if len(connected) == 0:
                current_peak = BinaryContinousAreas(self.graph)
                current_peak.add(single_area.node, single_area.start,
                                 single_area.end)
                if single_area.is_start and single_area.is_end:
                    peaks.append(current_peak)
                    #print("   Creating new, finishing it")
                else:
                    active_subgraphs.append(current_peak)
                    #print("   Creating new, appending to active")
            else:
                #print("   Adding to active")
                current_peak = connected[0]

                current_peak.add(single_area.node,
                             single_area.start,
                             single_area.end)

            # If more than 2, the first must be finished
            if len(active_subgraphs) > 5:
                peaks.append(active_subgraphs.popleft())
                #print("Finishing subgraph")

        peaks.extend(active_subgraphs)

        return peaks


    def build(self):
        self.single_areas = self._make_single_areas()
        return self._create_subgraphs()

    """
    def build(self):
        sorted_nodes = sorted(self.areas.keys())

        starts = OrderedDict()
        ends = OrderedDict()

        for node in sorted_nodes:
            starts_ends = self.areas[node]
            assert len(starts_ends) == 2 or len(starts_ends) == 4

            if self.areas._is_start_position(node, starts_ends[0]):
                starts[node] = starts_ends[:2]

            if self.areas._is_end_position(node, starts_ends[-1]):
                ends[node] = starts_ends[-2:]


        # Iterate sorted nodes. Count number of starts and number of ends
        in_peak = False
        prev_status = in_peak
        peaks = []
        current_peak = BinaryContinousAreas(self.graph)
        for node in sorted_nodes:
            start_ends = self.areas[node]
            assert len(start_ends) == 2 or len(starts_ends) == 4
            is_start = node in starts
            is_end = node in ends

            if is_start and is_end:
                current_peak.add(node, start_ends[0:2])
    """
