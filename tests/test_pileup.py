import unittest
from graph_peak_caller.pileup import Pileup
from examples import *


class TestPileup(unittest.TestCase):
    def test_create(self):
        pileup = Pileup(ob_graphs[0], pileup_intervals)
        pileup.create()
        self.assertEqual(pileup, true_pileup)

if __name__ == "__main__":
    unittest.main()
