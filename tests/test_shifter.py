import unittest
from graph_peak_caller.shifter import Shifter
from examples import *
from offsetbasedgraph import Interval, Block, Graph



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
    def setUp(self):
        self.nodes = {i: Block(1000) for i in range(10)}
        self.edges = {i: [i+1] for i in range(9)}
        self.graph = Graph(self.nodes, self.edges)
        self.interval = Interval(400, 600, [5])
        self.neg_interval = Interval(400, 600, [5], direction=-1)
        self.shifter = Shifter(self.graph, 1000)
        self.small_shifter = Shifter(self.graph, 300)

    def test_extend(self):
        areas = self.shifter.extend_interval_fast(self.interval)
        true_areas = {5: [400, 1000], 6: [0, 400]}
        self.assertEqual(areas, true_areas)

    def test_extend_small(self):
        areas = self.small_shifter.extend_interval_fast(self.interval)
        true_areas = {5: [400, 700]}
        self.assertEqual(areas, true_areas)

    def test_extend_neg_small(self):
        areas = self.small_shifter.extend_interval_fast(self.neg_interval)
        true_areas = {5: [300, 600]}
        self.assertEqual(areas, true_areas)

    def test_extend_neg(self):
        areas = self.shifter.extend_interval_fast(self.neg_interval)
        true_areas = {4: [600, 1000], 5: [0, 600]}
        self.assertEqual(areas, true_areas)

    def test_extend_both(self):
        areas = self.shifter.extend_interval_fast(self.interval, 0)
        true_areas = {4: [400, 1000], 5: [0, 1000], 6: [0, 400]}
        self.assertEqual(areas, true_areas)

    def _test_extend_interval(self):
        graph = ob_graphs[0]
        interval = pileup_intervals[0]
        d = 3
        shifter = Shifter(graph, pileup_intervals, d)
        areas = shifter.extend_interval(interval)
        true_areas = {0: [5, 10],
                      1: [0, 8]}
        self.assertEqual(areas, true_areas)

    def _test_extend_branch(self):
        shifter = Shifter(ob_graphs[0], [], 5)
        interval = Interval(0, 10, [0], graph=ob_graphs[0])
        areas = shifter.extend_interval(interval)
        true_areas = {0: [0, 10],
                      1: [0, 5],
                      2: [0, 5]}

        self.assertEqual(areas, true_areas)

    def _test_extend_disjoint(self):
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

    def _test_extend_disjoint_reverse(self):
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

    def _test_extend_both_sides(self):
        graph = ob_graphs[0].copy()
        graph.blocks[2]._length = 22
        shifter = Shifter(graph, [], 10)
        interval = Interval(5, 15, [1], graph=graph, direction=-1)
        areas = shifter.extend_interval(interval, 0)
        true_areas = {
            3: [0, 5],
            1: [0, 20],
            0: [5, 10]}
        print("#####", areas)
        self.assertEqual(areas, true_areas)


if __name__ == "__main__":
    unittest.main()
