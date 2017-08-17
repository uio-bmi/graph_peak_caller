import unittest
from graph_peak_caller.shifter import Shifter
from offsetbasedgraph import Interval
from examples import *


class TestShifter(unittest.TestCase):
    def _test_shift_interval(self):
        graph = ob_graphs[0]
        interval = pileup_intervals[0]
        d = 3
        shifter = Shifter(graph, pileup_intervals, d)
        areas = shifter.shift_interval(interval)
        true_areas = {0: [8, 10],
                      1: [0, 8]}
        self.assertEqual(areas, true_areas)

    def test_shift_branch(self):
        shifter = Shifter(ob_graphs[0], [], 5)
        interval = Interval(0, 10, [0], graph=ob_graphs[0])
        areas = shifter.shift_interval(interval)
        true_areas = {0: [5, 10],
                      1: [0, 5],
                      2: [0, 5]}
        self.assertEqual(areas, true_areas)

    def test_shift_disjoint(self):
        graph = ob_graphs[0].copy()
        graph.blocks[2]._length = 22
        shifter = Shifter(graph, [], 25)
        interval = Interval(0, 10, [0], graph=graph)
        areas = shifter.shift_interval(interval)
        true_areas = {1: [15, 20],
                      2: [15, 22],
                      3: [0, 5]}

        self.assertEqual(areas, true_areas)


if __name__ == "__main__":
    unittest.main()
