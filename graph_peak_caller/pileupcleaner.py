from itertools import chain
import numpy as np
import offsetbasedgraph as obg
from .extender import Areas
from collections import defaultdict
from graph_peak_caller.sparsepileup import SparsePileup

class IntervalWithinBlock(obg.Interval):
    def __init__(self, id, start, end, region_paths, graph):
        #assert len(region_paths) == 1
        super(IntervalWithinBlock, self).__init__(start, end, region_paths, graph)
        self.id = id
        self.deleted = False
        self.has_infinite_loop = False
        self.has_been_expanded = False
        self.direction = False

    def is_at_beginning_of_block(self):
        return self.start_position.offset == 0

    def is_at_end_of_block(self):
        return self.end_position.offset == \
               self.graph.blocks[self.end_position.region_path_id].length()

    def blocks_going_into(self):
        if not self.is_at_beginning_of_block():
            return []
        else:
            return self.graph.reverse_adj_list[-self.region_paths[0]]

    def blocks_going_out_from(self):
        if not self.is_at_end_of_block():
            return []
        else:
            end_rp = self.region_paths[-1]
            return self.graph.adj_list[end_rp] + \
                   list(np.abs(self.graph.reverse_adj_list[end_rp]))

    def merge_right(self, other):
        direction = False
        if other.region_paths[0] in self.graph.adj_list[self.region_paths[-1]]:
            direction = 1
        elif -self.region_paths[-1] in self.graph.adj_list[-other.region_paths[0]]:
            direction = -1
        else:
            raise Exception("Cannot merge")

        if self.direction and self.direction != direction:
            print("Not merging, conflicting direction")
            return False

        interval = IntervalWithinBlock(-1, self.start_position, other.end_position,
                                   self.region_paths + other.region_paths,
                                   self.graph)
        interval.direction = direction
        return interval

class PileupCleaner(object):


    def __init__(self, pileup):
        self.pileup  = pileup
        self.graph = pileup.graph
        self.valued_areas = self.pileup.find_valued_areas(True)
        self.non_valued_areas = self.pileup.find_valued_areas(False)

        self.intervals_at_start_of_block = defaultdict(list)
        self.intervals_at_end_of_block = defaultdict(list)

        self.intervals = []

    def remove_holes(self, threshold):
        holes_intervals = self.find_trivial_intervals_within_blocks(self.non_valued_areas)


    def filter_on_length_and_return_pileup(self, threshold):
        self.find_trivial_intervals_within_blocks(self.valued_areas)
        filtered_intervals = self.filter_on_length(threshold)
        return SparsePileup.from_intervals(self.graph, filtered_intervals)

    def filter_on_length(self, threshold):
        self.create_interval_indices()
        self._remove_short_intervals(threshold)

        while True:
            print("=== Mergin with next ==")
            if not self._merge_intervals_with_next_blocks():
                break

        # Filter
        filtered_intervals = []
        for interval in self.intervals:
            if interval.deleted:
                continue

            if interval.length() >= threshold or interval.has_infinite_loop:
                filtered_intervals.append(interval)

        return filtered_intervals


    def get_long_intervals(self, length_threshold):
        pass

    def find_trivial_intervals_within_blocks(self, areas):
        intervals = []
        id_counter = 0


        for node in areas:
            #print("Node %d" % node)
            #print(self.valued_areas[node])
            for i in range(0, len(areas[node]) // 2):
                start = int(areas[node][i*2])
                end = int(areas[node][i*2+1])
                #print("Start/end %d/%d" % (start, end))
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

    def _merge_single_interval_with_nexts(self, interval):
        # Find candidates to merge
        merge_with_intervals = []
        #for next_block in self.graph.adj_list[interval.region_paths[-1]]:
        for next_block in interval.blocks_going_out_from():
            print("     Checking block %d" % next_block)
            for next_interval in self.intervals_at_start_of_block[next_block]:
                print("    Merging with %s" % str(next_interval))

                if next_interval.region_paths[0] == interval.region_paths[0]:
                    print("Loop detected")
                    interval.has_infinite_loop = True
                    continue

                if not next_interval.deleted:
                    merge_with_intervals.append(next_interval)
                else:
                    print("    Deleted")


        if len(merge_with_intervals) == 0:
            return False

        n_new_intervals = 0
        interval.has_been_expanded = True  # Delete original,  make new
        for other_interval in merge_with_intervals:
            #other_interval.deleted = True

            new_interval = interval.merge_right(other_interval)
            if new_interval:
                new_interval.id = len(self.intervals)
                self.intervals.append(new_interval)
                n_new_intervals += 1

        if n_new_intervals:
            return True
        return False

    def _merge_intervals_with_next_blocks(self):
        n_merged = 0
        for interval in self.intervals:
            print("  Mering interval %s" % str(interval))
            if interval.deleted or interval.has_been_expanded:
                continue

            if not interval.is_at_end_of_block():
                continue

            merge_result = self._merge_single_interval_with_nexts(interval)
            if merge_result:
                n_merged += 1

        if n_merged > 0:
            return True

        return False

    def create_interval_indices(self):
        # Create indices from interval id to touching blocks
        for interval in self.intervals:
            if interval.is_at_end_of_block():
                self.intervals_at_end_of_block[interval.region_paths[-1]].append(interval)
            if interval.is_at_beginning_of_block():
                self.intervals_at_start_of_block[interval.region_paths[0]].append(interval)



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

