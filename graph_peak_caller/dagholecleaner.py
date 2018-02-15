from .densepileup import DensePileup
import numpy as np
from .extender import Extender
import offsetbasedgraph as obg
import logging


class DagHoleCleaner(object):
    """
    Efficient Hole cleaner for DensePileup on Directed Acyclic Graph
    """

    def __init__(self, dense_pileup, hole_size):
        assert isinstance(dense_pileup, DensePileup)
        self.pileup = dense_pileup
        self._graph = dense_pileup.graph
        self.hole_size = hole_size

    def run(self):
        positions = self.get_left_side_of_holes()
        self.expand_hole_sides_to_right(positions)
        return self.pileup

    def get_left_side_of_holes(self):
        """
        Find all internal left sides by checking diff == -1
         and ignoring beginning of nodes
        Add all left hole side that are at beginning of node by
        explicitly checking all node starts
        """
        logging.info("Getting left side of holes")
        node_lengths = [self.pileup.data.node_size(node)
                        for node in self.pileup.data._nodes[0:-1]]
        start_node_indices = np.cumsum(node_lengths)
        diffs = np.diff(self.pileup.data._values * 1)

        if len(start_node_indices) > 0:
            diffs[start_node_indices-1] = 0  # Ignore start of nodes here, we cannot be sure they are holes

        left_hole_indices = np.where(diffs == -1)[0] + 1
        positions = self.pileup.data.value_indexes_to_nodes_and_offsets(
            left_hole_indices)

        # Check all node starts
        logging.info("Adding node starts")
        is_hole = set(
            np.where(self.pileup.data._values == False)[0])
        logging.info("Step 1 done")
        is_hole_and_start = is_hole.intersection(set(start_node_indices))
        logging.info("Step 2 done")
        is_hole_and_start = sorted(list(is_hole_and_start))
        logging.info("Step 3 done")
        start_positions = self.pileup.data.value_indexes_to_nodes_and_offsets(
            is_hole_and_start)

        logging.info("Going through each start")
        for start_position in start_positions:
            node = start_position[0]
            is_start_of_hole = True
            if len(self._graph.reverse_adj_list[-node]) == 0:
                continue  # Start of grap, not hole

            for in_node in self._graph.reverse_adj_list[-node]:
                values = self.pileup.data.values(-in_node)
                if not values[-1]:
                    is_start_of_hole = False
                    break

            if is_start_of_hole:
                positions.append(start_position)

        return positions

    def expand_hole_sides_to_right(self, positions):

        logging.info("Number of holes to expand: %d" % len(positions))
        extender = Extender(self._graph, self.hole_size)

        i = 0
        for position in positions:
            if i % 500 == 0:
                logging.info("Expanding hole %d/%d" % (i, len(positions)))
            i += 1

            node = position[0]
            offset = position[1]
            interval = obg.DirectedInterval(
                int(offset), int(offset) + 1, [node])
            areas = extender.extend_interval(interval)
            self.pileup.set_area_to_value(areas, True)

        self.pileup.threshold(0.5)
