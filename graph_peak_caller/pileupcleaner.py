from itertools import chain
import numpy as np
import offsetbasedgraph as obg
from .extender import Areas
from collections import defaultdict
import logging

logging.basicConfig(level=logging.ERROR)


class IntervalWithinBlock(obg.Interval):
    def __init__(self, id, start, end, region_paths, graph):
        #assert len(region_paths) == 1
        super(IntervalWithinBlock, self).__init__(start, end, region_paths, graph)
        self.id = id
        self.deleted = False
        self.has_infinite_loop = False
        self.has_been_expanded = False
        self.direction = False
        self.is_maximally_expanded = False
        self.n_possible_expansions = 0
        self.n_expansions = 0

    def is_at_beginning_of_block(self):
        return self.start_position.offset == 0

    def is_at_end_of_block(self):
        return self.end_position.offset == \
               self.graph.blocks[self.end_position.region_path_id].length()

    def blocks_going_into(self):
        if not self.is_at_beginning_of_block():
            return []
        else:
            return list(np.abs(self.graph.reverse_adj_list[-self.region_paths[0]])) + \
                    list(np.abs(self.graph.adj_list[-self.region_paths[0]]))

    def blocks_going_out_from(self):
        if not self.is_at_end_of_block():
            return []
        else:
            end_rp = self.region_paths[-1]
            return self.graph.adj_list[end_rp] + \
                   list(np.abs(self.graph.reverse_adj_list[end_rp]))

    def merge(self, other):
        self.n_expansions += 1
        #other.n_expansions += 1

        # Check if we are merging on left or right side
        if other.region_paths[0] in self.blocks_going_out_from():
            return self.merge_right(other)
        elif self.region_paths[0] in other.blocks_going_out_from():
            return other.merge_right(self)
        else:
            raise Exception("Trying to merge intervals not connected")

    def merge_right(self, other):
        if other.region_paths[0] in self.graph.adj_list[self.region_paths[-1]]:
            direction = 1
        elif -self.region_paths[-1] in self.graph.adj_list[-other.region_paths[0]]:
            direction = -1
        else:
            raise Exception("Cannot merge, trying to merge with non-connecting interval")
        if (self.direction and self.direction != direction) or (other.direction and other.direction != direction):
            self.is_maximally_expanded = True
            return False

        interval = IntervalWithinBlock(-1, self.start_position, other.end_position,
                                   self.region_paths + other.region_paths,
                                   self.graph)
        interval.direction = direction
        return interval


class PileupCleaner(object):

    def __init__(self, pileup):
        self.pileup = pileup
        self.graph = pileup.graph
        self.valued_areas = self.pileup.find_valued_areas(True)
        logging.debug(self.valued_areas)
        self.non_valued_areas = self.pileup.find_valued_areas(False)
        self.intervals_at_start_of_block = defaultdict(list)
        self.intervals_at_end_of_block = defaultdict(list)

        self.intervals = []

    def get_small_holes(self, threshold):
        maximally_expanded_holes = self.find_maximally_expanded_holes(threshold)
        filtered = []
        for interval in maximally_expanded_holes:
            if interval.length() <= threshold:
                filtered.append(interval)

        return filtered

    def find_maximally_expanded_holes(self, max_size=False):
        """
        Algorithm:
         - Expand every hole maximally
         - In the end, only consider holes that have been maximally expanded
         - These are either intervals has_been_expanded = False, OR
           intervals that have been expanded fewer times than possible paths to expand
            (latter meaning that it has been maximally expanded in one direction,
            since it could not have been exdended further in that direction
        """

        self.intervals = self.find_trivial_intervals_within_blocks(self.non_valued_areas)
        self.create_interval_indices()

        self._set_number_of_possible_expansions_for_intervals()

        while True:
            if not self._merge_intervals_with_next_blocks(max_size):
                break

        maximally_expanded = []
        for interval in self.intervals:
            if not interval.has_been_expanded \
                    or interval.n_expansions < interval.n_possible_expansions:
                maximally_expanded.append(interval)

        return maximally_expanded

    def filter_on_length(self, threshold):
        logging.debug("filter_on_length: %s", threshold)
        self.find_trivial_intervals_within_blocks(self.valued_areas)
        self.create_interval_indices()
        self._remove_short_intervals(threshold)

        while True:
            if not self._merge_intervals_with_next_blocks(threshold):
                break

        # Filter
        filtered_intervals = []
        for interval in self.intervals:
            if interval.deleted:
                continue

            if interval.length() >= threshold or interval.has_infinite_loop:
                filtered_intervals.append(interval)

        return filtered_intervals

    def _set_number_of_possible_expansions_for_intervals(self):
        for interval in self.intervals:
            interval.n_possible_expansions = 0
            if interval.is_at_beginning_of_block():
                interval.n_possible_expansions += len(interval.blocks_going_into())

            if interval.is_at_end_of_block():
                interval.n_possible_expansions += len(interval.blocks_going_out_from())


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
        if region_path in interval.region_paths[:-1]:
            return False
        return region_path == interval.region_paths[-1] and \
            interval.start_position.offset == self.graph.node_size(region_path)

    def _merge_single_interval_with_nexts(self, interval):
        # Find candidates to merge
        merge_with_intervals = []
        for next_block in interval.blocks_going_out_from():
            for next_interval in self.intervals_at_start_of_block[next_block]:

                if self.has_interval_region_path_start(
                        interval,
                        next_interval.region_paths[0]):
                    interval.has_infinite_loop = True
                    continue

                if not next_interval.deleted:
                    merge_with_intervals.append(next_interval)

        for next_block in interval.blocks_going_into():
            for next_interval in self.intervals_at_end_of_block[next_block]:
                if self.has_interval_region_path_end(
                        interval, next_interval.region_paths[-1]):
                    interval.has_infinite_loop = True
                    continue

                if not next_interval.deleted:
                    merge_with_intervals.append(next_interval)

        if len(merge_with_intervals) == 0:
            return False

        n_new_intervals = 0
        interval.has_been_expanded = True  # Delete original,  make new
        for other_interval in merge_with_intervals:
            #other_interval.deleted = True

            new_interval = interval.merge(other_interval)
            if new_interval:
                new_interval.id = len(self.intervals)
                self.intervals.append(new_interval)
                n_new_intervals += 1

        if n_new_intervals > 0:
            return True
        return False

    def _merge_intervals_with_next_blocks(self, max_length=False):
        n_merged = 0
        i = 0
        length = len(self.intervals)
        for interval in self.intervals[0: length]:
            if i % 1000 == 0:
                print("Merging interval %d / %d"  % (i, length)) 
            i += 1

            if interval.deleted or interval.has_been_expanded:
                continue

            if not interval.is_at_end_of_block() and not interval.is_at_beginning_of_block():
                continue

            if max_length and interval.length() > max_length:
                interval.has_been_expanded = True
                continue

            merge_result = self._merge_single_interval_with_nexts(interval)
            if merge_result:
                n_merged += 1

        if n_merged > 0:
            return True

        return False

    def create_interval_indices(self):
        # Create indices from interval id to touching blocks
        logging.debug("create_interval_indices")
        for interval in self.intervals:
            logging.debug("Testing: %s", interval)
            if interval.is_at_end_of_block():
                logging.debug("End of Block: %s", interval)
                self.intervals_at_end_of_block[interval.region_paths[-1]].append(interval)
            if interval.is_at_beginning_of_block():
                logging.debug("Beggining of Block: %s", interval)
                self.intervals_at_start_of_block[interval.region_paths[0]].append(interval)


    # Delete
    def _get_inverse_valued_areas(self):
        inverse = {}
        for node_id in self.valued_areas:
            inverse_areas = []
            inverse_areas.append(0)
            inverse_areas.append(self.valued_areas)
            inverse_areas.append(self.graph.blocks[node_id].length())

            if inverse_areas[0] == inverse_areas[1]:
                inverse_areas = inverse_areas[2:]

            if inverse_areas[-2] == inverse_areas[-1]:
                inverse_areas= inverse_areas[0:-2]

            inverse[node_id] = inverse_areas

        return inverse

