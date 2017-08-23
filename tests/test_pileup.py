import unittest
from graph_peak_caller.pileup import Pileup
from examples import *


class TestPileup(unittest.TestCase):
    def test_create(self):
        pileup = Pileup(ob_graphs[0])
        pileup.add_intervals(pileup_intervals)
        self.assertEqual(pileup, true_pileup)

    def test_to_file_from_file(self):
        pileup =  Pileup(ob_graphs[0])
        pileup.add_intervals(pileup_intervals)
        pileup.to_bed_graph("test.pileup")

        pileup_from_file = Pileup.from_bed_graph(ob_graphs[0], "test.pileup")

        self.assertTrue(pileup_from_file == pileup)

    def test_fill_small_holes(self):
        pass

if __name__ == "__main__":
    unittest.main()
