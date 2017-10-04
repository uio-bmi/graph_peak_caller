from itertools import chain
import numpy as np
import offsetbasedgraph as obg
from .extender import Areas
from collections import defaultdict
import logging
from .pileupcleaner import IntervalWithinBlock

logging.basicConfig(level=logging.ERROR)

class PileupCleaner2(object):

    def __init__(self, pileup):
        self.pileup = pileup
        self.graph = pileup.graph
        self.valued_areas = self.pileup.find_valued_areas(True)
        logging.debug(self.valued_areas)
        self.non_valued_areas = self.pileup.find_valued_areas(False)
        self.intervals_at_start_of_block = defaultdict(list)
        self.intervals_at_end_of_block = defaultdict(list)
        self._merged_intervals_counter = 0
        self.finding_holes = False

        self.intervals = []

    def get_small_holes(self, threshold):
        self.intervals = self.find_trivial_intervals_within_blocks(self.non_valued_areas)
        self.create_interval_indices()
        self.merge_intervals_iteratively()

        return self.get_intervals_from_indices_above_or_below_length(threshold, above=False)

    def filter_on_length(self, threshold):

        self.intervals = self.find_trivial_intervals_within_blocks(self.valued_areas)
        self.create_interval_indices()
        self.merge_intervals_iteratively()

        return self.get_intervals_from_indices_above_or_below_length(threshold, above=True)


    def get_intervals_from_indices_above_or_below_length(self, threshold, above=True):
        filtered = []
        # Include all long enough trivial intervals
        for interval in self.intervals:
            if above and interval.length() >= threshold:
                filtered.append(interval)
            if not above and interval.length() < threshold:
                filtered.append(interval)

        for node, intervals in self.intervals_at_end_of_block.items():
            for interval in intervals:
                if interval not in filtered:
                    if above and (interval.length() >= threshold or interval.has_infinite_loop):
                        filtered.append(interval)
                    if not above and interval.length() < threshold:
                        filtered.append(interval)

        for node, intervals in self.intervals_at_start_of_block.items():
            for interval in intervals:
                if interval not in filtered:
                    if above and (interval.length() >= threshold or interval.has_infinite_loop):
                        filtered.append(interval)
                    if not above and interval.length() < threshold:
                        filtered.append(interval)

        logging.debug("== Filtered ==")
        for interval in filtered:
            logging.debug(" %s" % str(interval))
        return filtered


    def merge_intervals_iteratively(self):

        while self._merge_intervals(1):
            logging.debug("== merged forward ==")

    def _merge_intervals(self, direction=1):
        print("Mering intervals")
        logging.debug(" MERGING INTERVALS (%d)" % direction)
        n_merged = 0
        nodes = list(self.intervals_at_end_of_block)
        i = 0
        for node_id in nodes:
            intervals = self.intervals_at_end_of_block[node_id]
            if i % 5000 == 0 or i < 200:
                n_intervals = sum([len(i) for i in self.intervals_at_end_of_block.values()])
                print("Merging node %d / %d (%d intervals, total %d)" % (i, len(nodes), len(intervals), n_intervals))
            i += 1
            logging.debug("   Merging for node %d " % (node_id))
            #print(intervals)
            for interval in intervals:
                if interval.has_been_expanded:
                    logging.debug("   has been expanded, skipping")
                    continue
                #logging.debug("        Merging interval %s" % interval)
                did_merge = self._merge_interval_right(interval, direction=direction)
                #print("intervals after merge")
                #print(self.intervals_at_end_of_block[node_id])
                self._remove_interval_from_end_indices(interval)
                #if self._n_intervals_into_node(interval.region_paths[0]) == 0:
                #           self._remove_interval_from_start_indices(interval)
                if did_merge:
                    n_merged += 1

                interval.has_been_expanded = True

        logging.debug("           n merged: %d" % n_merged)
        if n_merged > 0:
            return True
        return False

    def _merge_interval_right(self, interval, direction):
        outgoing_edges = interval.blocks_going_out_from(limit_to_direction=direction)
        #logging.debug("       Outgoing edges: %s " % str(outgoing_edges))
        n_merged = 0
        for node in outgoing_edges:
            #logging.debug("           Checking node %d" % node)
            touching_intervals = self.intervals_at_start_of_block[node]
            #logging.debug("            Found touching: %s" % str(touching_intervals))
            for touching_interval in touching_intervals:
                #logging.debug("            Merging with %s" % str(touching_interval))
                if self.is_forward_merging_into_loop(interval, touching_interval):
                    logging.debug("              LOOP DETECTED")
                    interval.has_infinite_loop = True
                    continue

                merged_interval = interval.merge(touching_interval)
                self._add_indices_for_interval(merged_interval)
                # If the one we merge into has only one edge in, we can remove it from the indices
                n_edges_in = len(self.graph.reverse_adj_list[-touching_interval.region_paths[0]])
                n_intervals_in = self._n_intervals_into_node(touching_interval.region_paths[0])
                logging.debug("                 n edges in: %d" % n_edges_in)
                logging.debug("                 n intervals in: %d" % n_intervals_in)
                if n_intervals_in == 1:  # 1 interval in meaning no more intervals in than the one currently being merged
                    logging.debug("!!!!            No more intervals possibly in, removing from indices")
                    self._remove_interval_from_start_indices(touching_interval)
                    self._remove_interval_from_end_indices(touching_interval)

                n_merged += 1

        if n_merged > 0:
            return True
        logging.debug("              No merging")
        return False

    def _n_intervals_into_node(self, node_id):
        n = 0
        nodes_in = self.graph.reverse_adj_list[-node_id]
        for node in nodes_in:
            n += len(self.intervals_at_end_of_block[-node])
        return n

    def _merge_interval_left(self, interval, direction):
        # Find all incoming edges, merge with intervals touching, update indices
        incoming_nodes = interval.blocks_going_into(limit_to_direction=direction)
        for node in incoming_nodes:
            touching_intervals = self.intervals_at_end_of_block[node]
            for touching_interval in touching_intervals:
                merged_interval = interval.merge(touching_interval)
                self._add_indices_for_interval(merged_interval)

        # Remove interval from start index. We do not want to expand it left again
        self._remove_interval_from_start_indices(interval)

    def _remove_interval_from_start_indices(self, interval):
        #logging.debug("            Removing from start: %s" % str(interval))
        start_rp = interval.region_paths[0]

        if interval in self.intervals_at_start_of_block[start_rp]:
            self.intervals_at_start_of_block[start_rp].remove(interval)

        #if interval in self.intervals_at_end_of_block[-start_rp]:
        #    self.intervals_at_end_of_block[-start_rp].remove(interval)

    def _remove_interval_from_end_indices(self, interval):
        #logging.debug("            Removing from end %s" % str(interval))
        end_rp = interval.region_paths[-1]
        if interval in self.intervals_at_end_of_block[end_rp]:
            self.intervals_at_end_of_block[end_rp].remove(interval)
            #print("Removed")
        #if interval in self.intervals_at_start_of_block[-end_rp]:
        #    self.intervals_at_start_of_block[-end_rp].remove(interval)
            #print("Removed")
        #print(self.intervals_at_end_of_block[end_rp])

    def _add_indices_for_interval(self, interval):
        start_rp = interval.region_paths[0]
        end_rp = interval.region_paths[-1]
        self.intervals_at_start_of_block[start_rp].append(interval)
        #self.intervals_at_end_of_block[-start_rp].append(interval)
        #self.intervals_at_start_of_block[-end_rp].append(interval)
        self.intervals_at_end_of_block[end_rp].append(interval)

    def _set_number_of_possible_expansions_for_intervals(self):
        for interval in self.intervals:
            interval.n_possible_expansions = 0
            if interval.is_at_beginning_of_block():
                interval.n_possible_expansions += len(interval.blocks_going_into())

            if interval.is_at_end_of_block():
                interval.n_possible_expansions += len(interval.blocks_going_out_from())


    def is_forward_merging_into_loop(self, interval, next_interval):
        if interval.contains(next_interval):
            return True

    def get_long_intervals(self, length_threshold):
        pass

    def find_trivial_intervals_within_blocks(self, areas):
        intervals = []
        id_counter = 0

        for node in areas:
            for i in range(0, len(areas[node]) // 2):
                start = int(areas[node][i*2])
                end = int(areas[node][i*2+1])
                interval = IntervalWithinBlock(id_counter, start, end, [node], self.graph)
                intervals.append(interval)
                interval.id = id_counter
                id_counter += 1

                node_size = self.graph.node_size(node)
                interval_reverse = IntervalWithinBlock(id_counter, node_size - end, node_size - start, [-node], self.graph)
                intervals.append(interval_reverse)
                id_counter += 1

        self.intervals = intervals
        return intervals

    def _remove_short_intervals(self, threshold):
        for interval in self.intervals:
            if interval.length() < threshold and  \
                    not interval.is_at_end_of_block() and \
                    not interval.is_at_beginning_of_block():
                interval.removed = True

    def has_interval_region_path_start(self, interval, region_path):
        if region_path in interval.region_paths[1:]:
            return False
        return region_path == interval.region_paths[0] and \
            interval.start_position.offset == 0

    def has_interval_region_path_end(self, interval, region_path):
        if region_path in interval.region_paths[:-1] and region_path != interval.region_paths[0]:
            return False
        return region_path == interval.region_paths[-1] and \
            interval.start_position.offset == self.graph.node_size(region_path)

    def _merge_single_interval_with_nexts(self, interval):
        pass

    def create_interval_indices(self):
        print(" === Create interval indices  === ")
        # Create indices from interval id to touching blocks
        logging.debug("create_interval_indices")
        for interval in self.intervals:
            #logging.debug("Testing: %s", interval)
            if interval.is_at_end_of_block():
                #logging.debug("End of Block: %s", interval)
                if interval not in self.intervals_at_end_of_block[interval.region_paths[-1]]:
                    self.intervals_at_end_of_block[interval.region_paths[-1]].append(interval)
                    #self.intervals_at_start_of_block[-interval.region_paths[-1]].append(interval)
            if interval.is_at_beginning_of_block():
                #logging.debug("Beginning of Block: %s", interval)
                if interval not in self.intervals_at_start_of_block[interval.region_paths[0]]:
                    self.intervals_at_start_of_block[interval.region_paths[0]].append(interval)
                    #self.intervals_at_end_of_block[-interval.region_paths[0]].append(interval)
