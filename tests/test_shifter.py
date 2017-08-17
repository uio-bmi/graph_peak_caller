import unittest
from graph_peak_caller.shifter import Shifter
from examples import *


class TestShifter(unittest.TestCase):
    def test_shift_interval(self):
        graph = ob_graphs[0]
        interval = pileup_intervals[0]
        d = 3
        shifter = Shifter(graph, pileup_intervals, d)
        areas = shifter.shift_interval(interval)
        true_areas = {0: [8, 10],
                      1: [0, 8]}
        self.assertEqual(areas, true_areas)

if __name__ == "__main__":
    unittest.main()
