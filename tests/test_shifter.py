import unittest
from graph_peak_caller.shifter import Shifter
from offsetbasedgraph import Interval
from examples import *


# class TestShifter(unittest.TestCase):
class TestShifter:
    def test_shift_interval(self):
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

    def test_shift_disjoint_reverse(self):
        graph = ob_graphs[0].copy()
        graph.blocks[2]._length = 22
        shifter = Shifter(graph, [], -25)
        interval = Interval(0, 10, [3], graph=graph)
        areas = shifter.shift_interval(interval)
        true_areas = {1: [0, 5],
                      2: [0, 7],
                      0: [5, 10]}
        print("#####", areas)
        self.assertEqual(areas, true_areas)


class TestExtend(unittest.TestCase):
    def test_extend_interval(self):
        graph = ob_graphs[0]
        interval = pileup_intervals[0]
        d = 3
        shifter = Shifter(graph, pileup_intervals, d)
        areas = shifter.extend_interval(interval)
        true_areas = {0: [5, 10],
                      1: [0, 8]}
        self.assertEqual(areas, true_areas)

    def test_extend_branch(self):
        shifter = Shifter(ob_graphs[0], [], 5)
        interval = Interval(0, 10, [0], graph=ob_graphs[0])
        areas = shifter.extend_interval(interval)
        true_areas = {0: [0, 10],
                      1: [0, 5],
                      2: [0, 5]}

        self.assertEqual(areas, true_areas)

    def test_extend_disjoint(self):
        graph = ob_graphs[0].copy()
        graph.blocks[2]._length = 22
        shifter = Shifter(graph, [], 25)
        interval = Interval(0, 10, [0], graph=graph)
        areas = shifter.extend_interval(interval)
        true_areas = {
            0: [0, 10],
            1: [0, 20],
            2: [0, 22],
            3: [0, 5]}

        self.assertEqual(areas, true_areas)

    def test_extend_disjoint_reverse(self):
        graph = ob_graphs[0].copy()
        graph.blocks[2]._length = 22
        shifter = Shifter(graph, [], 25)
        interval = Interval(0, 10, [3], graph=graph, direction=-1)
        areas = shifter.extend_interval(interval)
        true_areas = {
            3: [0, 10],
            1: [0, 20],
            2: [0, 22],
            0: [5, 10]}
        print("#####", areas)
        self.assertEqual(areas, true_areas)


if __name__ == "__main__":
    unittest.main()
