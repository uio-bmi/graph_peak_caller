import unittest
from graph_peak_caller.pileup import Pileup, SmallIntervals

from examples import *


class TestSmallIntervals(unittest.TestCase):
    def setUp(self):
        graph = ob_graphs[0]
        end_interval_dict = {0: (8, 10)}
        whole_intervals = [1, 2]
        interval_dict = {3: [(0, 3)]}
        self.small_intervals = SmallIntervals(interval_dict, end_interval_dict, whole_intervals, graph, 26)

    def test_get_intervals(self):
        intervals = self.small_intervals._get_intervals(0, [8])
        true_intervals = [[8, 0, 1, 3, 3],
                          [8, 0, 2, 3, 3]]
        self.assertEqual(intervals, true_intervals)

if __name__ == "__main__":
    unittest.main()
