from graph_peak_caller.densepileup import DensePileup
import numpy as np

class DagHoleCleaner(object):
    """
    Efficient Hole cleaner for DensePileup on Directed Acyclic Graph
    """

    def __init__(self, dense_pileup, hole_size):
        assert isinstance(dense_pileup, DensePileup)
        self.pileup = dense_pileup
        self._graph = dense_pileup.graph

    def run(self):
        positions = self.find_left_side_of_holes()
        self.expand_hole_sides_to_right(positions)
        self.expand_back()

    def find_left_side_of_holes(self):
        diffs = np.diff(self.pileup.data._values * 1)
        left_hole_indices = np.where(diffs == -1) + 1
        positions = self.pileup.data.value_indexes_to_nodes_and_offsets(left_hole_indices)
        return positions

    def expand_hole_sides_to_right(self, positions):


    def expand_back(self):
        pass

