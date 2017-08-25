import unittest
from graph_peak_caller.pileup import Pileup
from offsetbasedgraph import Graph, Block, Interval
from graph_peak_caller import bdgcmp
import numpy as np
from examples import *



class TestPileup(unittest.TestCase):
    def test_simple_bdgcmp(self):
        graph = one_block_graph
        pileup1 = pileup1_one_block
        pileup2 = pileup2_one_block

        max_pileup = bdgcmp.create_background_pileup_as_max_from_pileups(graph, iter([pileup1, pileup2]), 0)

        correct_count_array = np.array([1, 2, 2, 2, 2, 1, 1, 1, 1, 1])
        closeness = np.isclose(max_pileup.get_count_arrays()[1], correct_count_array)
        print(closeness)
        print(max_pileup.get_count_arrays()[1])
        self.assertTrue(np.all(np.isclose(max_pileup.get_count_arrays()[1], correct_count_array)))



if __name__ == "__main__":
    unittest.main()
