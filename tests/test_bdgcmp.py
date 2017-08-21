import unittest
from graph_peak_caller.pileup import Pileup
from offsetbasedgraph import Graph, Block, Interval
from bdgcmp import create_background_pileup_as_max_from_pileups
import numpy as np
from examples import *



class TestPileup(unittest.TestCase):
    def test_simple_bdgcmp(self):
        graph = one_block_graph
        pileup1 = pileup1_one_block
        pileup2 = pileup2_one_block

        max_pileup = create_background_pileup_as_max_from_pileups([pileup1, pileup2])

        correct_count_array = np.array([2, 2, 2, 2, 2, 1, 1, 1, 1])
        self.assertTrue(max_pileup.count_arrays[1] == correct_count_array)



if __name__ == "__main__":
    unittest.main()
