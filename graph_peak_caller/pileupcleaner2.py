from itertools import chain
import numpy as np
import offsetbasedgraph as obg
from .extender import Areas
from collections import defaultdict
import logging
from .pileupcleaner import IntervalWithinBlock, IndexedList

logging.basicConfig(level=logging.ERROR)

class PileupCleaner2(object):

    def __init__(self, pileup):
        self.pileup = pileup
        self.graph = pileup.graph
        self.valued_areas = self.pileup.find_valued_areas(True)
        #logging.debug(self.valued_areas)
        self.non_valued_areas = self.pileup.find_valued_areas(False)
        self.intervals_at_start_of_block = defaultdict(IndexedList)
        self.intervals_at_end_of_block = defaultdict(IndexedList)
        self._merged_intervals_counter = 0
        self.finding_holes = False
        self._n_loops_detected = 0

        self.intervals = []
        self.individual_intervals = []

    def get_small_holes(self, threshold):
        self.intervals = self.find_trivial_intervals_within_blocks(self.non_valued_areas)
        self.create_interval_indices()
        self.merge_intervals_iteratively()

        filtered = self.get_intervals_from_indices_above_or_below_length(threshold, above=False)
        intervals = PileupCleaner2._intervals_to_trivial_intervals(filtered)
        intervals = PileupCleaner2._remove_interval_duplicates(intervals)
        return intervals

    def filter_on_length(self, threshold, return_single_rp_intervals=False):

        self.intervals = self.find_trivial_intervals_within_blocks(self.valued_areas)
        self.create_interval_indices()
        self.merge_intervals_iteratively()

        filtered = self.get_intervals_from_indices_above_or_below_length(threshold, above=True)
        if not return_single_rp_intervals:
            return filtered

        intervals = PileupCleaner2._intervals_to_trivial_intervals(filtered)
        intervals = PileupCleaner2._remove_interval_duplicates(intervals)
        return intervals

    @staticmethod
    def _remove_interval_duplicates(intervals):
        print("Remove duplicates")
        index = {}

        out = []
        i = 0
        for interval in intervals:

            if i % 5000 == 0:
                print("Removing duplicates %d / %d" % (i, len(intervals)))
            i += 1
            if interval.length() == 0:
                continue
            hash = interval.hash(ignore_direction=True)
            if hash not in index:
                out.append(interval)
            index[hash] = True

        return out


    def get_intervals_from_indices_above_or_below_length(self, threshold, above=True):
        print("get_intervals_from_indices_above_or_below_length")
        filtered = []

        candidates = self.individual_intervals + \
                        self.intervals

        for interval in candidates:
            if above and interval.length() >= threshold:
                filtered.append(interval)
            if not above and interval.length() < threshold:
                filtered.append(interval)

        for node, intervals in self.intervals_at_end_of_block.items():
            for interval in intervals.get_elements():
                if above and (interval.length() >= threshold or interval.has_infinite_loop):
                    filtered.append(interval)
                if not above and interval.length() < threshold:
                    filtered.append(interval)

        for node, intervals in self.intervals_at_start_of_block.items():
            for interval in intervals.get_elements():
                if above and (interval.length() >= threshold or interval.has_infinite_loop):
                    filtered.append(interval)
                if not above and interval.length() < threshold:
                    filtered.append(interval)

        return PileupCleaner2._remove_interval_duplicates(filtered)

    @staticmethod
    def _intervals_to_trivial_intervals(intervals):
        print("Convert to trivial intervals")
        # Convert to only trivial intervals
        trivial_intervals = []
        i = 0
        for interval in intervals:
            if i % 5000 == 0:
                print("Converting  %d / %d" % (i, len(intervals)))
            i += 1

            graph = interval.graph
            rps = interval.region_paths
            if len(rps) == 1:
                trivial_intervals.append(interval)
            else:
                start_interval = obg.Interval(interval.start_position.offset, graph.node_size(rps[0]), [rps[0]], graph)
                end_interval = obg.Interval(0, interval.end_position.offset, [rps[-1]], graph)
                trivial_intervals.append(start_interval)
                trivial_intervals.append(end_interval)

                for rp in rps[1:-1]:
                    trivial_intervals.append(obg.Interval(0, graph.node_size(rp), [rp], graph))

        # Invert inverse
        for i, interval in enumerate(trivial_intervals):
            rp = interval.region_paths[0]
            if rp < 0:
                rp_size = graph.node_size(rp)
                new_start = rp_size - interval.end_position.offset
                new_end = rp_size - interval.start_position.offset
                new_interval = obg.Interval(new_start, new_end, [-rp], interval.graph)
                trivial_intervals[i] = new_interval

        return trivial_intervals

    def merge_intervals_iteratively(self):

        while self._merge_intervals(1):
            #logging.debug("== merged forward ==")
            continue

    def _merge_intervals(self, direction=1):
        print("Mering intervals")
        #logging.debug(" MERGING INTERVALS (%d)" % direction)
        n_merged = 0
        nodes = list(self.intervals_at_end_of_block)
        i = 0
        for node_id in nodes:
            intervals = self.intervals_at_end_of_block[node_id]
            if i % 20 == 0 or i < 30:
                n_intervals = sum([len(i) for i in self.intervals_at_end_of_block.values()])
                print("Merging node %d / %d (%d intervals, total %d, %d loops)" % \
                      (i, len(nodes), len(intervals), n_intervals, self._n_loops_detected))
            i += 1
            #logging.debug("   Merging for node %d " % (node_id))
            for interval in list(intervals.get_elements()).copy():
                if interval.has_been_expanded:
                    #logging.debug("   has been expanded, skipping")
                    continue
                #logging.debug("      Merging interval %s" % interval)
                did_merge = self._merge_interval_right(interval, direction=direction)
                if did_merge:
                    self._remove_interval_from_start_indices(interval)
                    self._remove_interval_from_end_indices(interval)
                    n_merged += 1

                interval.has_been_expanded = True
            #if i > 200:
            #    break
        #logging.debug("           n merged: %d" % n_merged)
        print("%d merged" % n_merged)
        if n_merged > 0:
            return True
        return False

    def _merge_interval_right(self, interval, direction):
        outgoing_edges = interval.blocks_going_out_from(limit_to_direction=direction)
        #logging.debug("       Outgoing edges: %s " % str(outgoing_edges))
        n_merged = 0
        n_intervals_checked = 0
        for node in outgoing_edges:
            #logging.debug("           Checking node %d" % node)
            touching_intervals = self.intervals_at_start_of_block[node]
            #logging.debug("            Found touching: %s" % str(touching_intervals))
            for touching_interval in list(touching_intervals.get_elements()).copy():
                n_intervals_checked += 1
                if n_intervals_checked > 100:
                    print(" Warning: %d intervals checked" % n_intervals_checked)

                #logging.debug("            Merging with %s" % str(touching_interval))
                if interval == touching_interval or (interval.contains_position(touching_interval.end_position) and \
                        (self.is_forward_merging_into_loop(interval, touching_interval) \
                        or interval.contains_loop()) \
                        or touching_interval.contains_loop()):
                    #logging.debug("              LOOP DETECTED")
                    interval.has_infinite_loop = True
                    self._n_loops_detected += 1
                    continue

                merged_interval = interval.merge(touching_interval)
                self._add_indices_for_interval(merged_interval)
                # If the one we merge into has only one edge in, we can remove it from the indices
                n_edges_in = len(self.graph.reverse_adj_list[-touching_interval.region_paths[0]])
                n_intervals_in = self._n_intervals_into_node(touching_interval.region_paths[0], [touching_interval, interval])
                if n_intervals_in == 0:
                    #logging.debug("!!!!            No more intervals possibly in, removing from indices")
                    self._remove_interval_from_start_indices(touching_interval)
                    self._remove_interval_from_end_indices(touching_interval)

                n_merged += 1

        if n_merged > 0:
            return True
        #logging.debug("              No merging")
        return False

    def _n_intervals_into_node(self, node_id, ignore_intervals=[]):
        n = 0
        nodes_in = self.graph.reverse_adj_list[-node_id]
        for node in nodes_in:
            for interval_in in self.intervals_at_end_of_block[-node].get_elements():
                if not interval_in in ignore_intervals:
                    n += 1
        return n


    def _remove_interval_from_start_indices(self, interval):
        start_rp = interval.region_paths[0]

        if interval in self.intervals_at_start_of_block[start_rp]:
            self.intervals_at_start_of_block[start_rp].remove(interval)

    def _remove_interval_from_end_indices(self, interval):
        end_rp = interval.region_paths[-1]
        if interval in self.intervals_at_end_of_block[end_rp]:
            self.intervals_at_end_of_block[end_rp].remove(interval)

    def _add_indices_for_interval(self, interval):
        start_rp = interval.region_paths[0]
        end_rp = interval.region_paths[-1]
        added = False
        #logging.debug("                 Adding indices for %s" % str(interval))
        if interval.start_position.offset == 0:
            self.intervals_at_start_of_block[start_rp].append(interval)
            added = True

        if interval.end_position.offset == self.graph.node_size(end_rp):
            self.intervals_at_end_of_block[end_rp].append(interval)

            added = True

        if not added:
            print("Did not add %s" % str(interval))
            self.individual_intervals.append(interval)

    def _set_number_of_possible_expansions_for_intervals(self):
        for interval in self.intervals:
            interval.n_possible_expansions = 0
            if interval.is_at_beginning_of_block():
                interval.n_possible_expansions += len(interval.blocks_going_into())

            if interval.is_at_end_of_block():
                interval.n_possible_expansions += len(interval.blocks_going_out_from())


    def is_forward_merging_into_loop(self, interval, next_interval):
        rps = interval.region_paths
        rps_next = next_interval.region_paths

        end_rp = rps_next[-1]
        end_offset = next_interval.end_position.offset

        for rp in rps_next:
            if rp in rps[1:]:
                return True
            if rp == rps[0] and (rp != rps_next[-1] or end_offset > interval.start_position.offset):
                return True

        return False

        # Loop if next interval ends inside current interval
        if not end_rp in rps:
            return False



        if end_rp == rps[0]:
            if end_offset > interval.start_position.offset:
                return True
            else:
                return False
        else:
            return True

        ###



        if len(rps) >= 3 and len(rps_next) == 1:
            if rps_next[0] in rps[1:-1]:
                return True

        if len(rps) > 1:
            if rps_next[0] in rps[1:]:
                return True
            if interval.start_position.offset == 0:
                if rps_next[0]  in rps:
                    return True

        if len(rps) == 1 and len(rps_next) > 1 and rps_next[0] == rps[0]:
            return True

        if next_interval.start_position.offset == 0 \
                and (rps_next[0] in rps[1:] or \
                interval.start_position.offset == 0 and rps_next[0] == rps[0]):
            return True

        if len(rps) > 1:
            if rps_next[-1] in rps[1:]:
                return True

        if len(rps) == 1 and len(rps_next) == 1 and rps[0] == rps_next[0]:
            if next_interval.end_position.offset > interval.start_position.offset:
                return True

        if interval.contains(next_interval):
            return True

        return False

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
        #logging.debug("create_interval_indices")
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

        self.sanitize_indices()

    def sanitize_indices(self):
        # Remove intervals that are within other intervals
        for index in [self.intervals_at_start_of_block, self.intervals_at_end_of_block]:
            for node in index:
                intervals = index[node]
                for interval in intervals.get_elements():
                    for other_interval in intervals.get_elements():
                        if other_interval == interval:
                            continue
                        if other_interval.contains(interval):
                            intervals.remove(interval)

