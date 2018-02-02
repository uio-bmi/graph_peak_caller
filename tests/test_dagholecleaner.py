from graph_peak_caller.densepileup import DensePileupData, DensePileup
from offsetbasedgraph import GraphWithReversals as Graph, Block, DirectedInterval as Interval
import unittest
import numpy as np
from graph_peak_caller.dagholecleaner import DagHoleCleaner

graph = Graph({i: Block(10) for i in range(1, 4)},
                                   {1: [2],
                                    2: [3]})

split_graph = Graph(
    {
        1: Block(10),
        2: Block(10),
        3: Block(10),
        4: Block(10),
    },
    {
        1: [2, 3],
        2: [4],
        3: [4]
    }
)

"""
class TestDagHoleCleanerGetLeftSideOfHoles(unittest.TestCase):

    def test_simple(self):
        pileup = DensePileup.from_intervals(graph,
                                            [Interval(0, 3, [1])]
                                            )

        cleaner = DagHoleCleaner(pileup, 3)
        left_holes = cleaner.get_left_side_of_holes()
        self.assertEqual(left_holes, [(1, 3)])

    def test_simple2(self):
        pileup = DensePileup.from_intervals(graph,
                                            [Interval(0, 3, [1])]
                                            )

        cleaner = DagHoleCleaner(pileup, 3)
        left_holes = cleaner.get_left_side_of_holes()
        self.assertEqual(left_holes, [(1, 3)])

    def test_special_case(self):
        pileup = DensePileup.from_intervals(graph,
                                            [Interval(0, 3, [1]),
                                             Interval(5, 10, [1])]
                                            )

        cleaner = DagHoleCleaner(pileup, 3)
        left_holes = cleaner.get_left_side_of_holes()
        self.assertEqual(left_holes, [(1, 3), (2, 0)])

    def test_special_case2(self):
        pileup = DensePileup.from_intervals(graph,
                                            [Interval(0, 3, [1]),
                                             Interval(5, 10, [2])]
                                            )

        cleaner = DagHoleCleaner(pileup, 3)
        left_holes = cleaner.get_left_side_of_holes()
        self.assertEqual(left_holes, [(1, 3), (3, 0)])


class TestDagHoleCleaner(unittest.TestCase):
    def test_single_peak(self):
        pileup = DensePileup.from_intervals(graph,
                                            [Interval(0, 3, [1])]
                                            )
        pileup.threshold(0.5)

        cleaner = DagHoleCleaner(pileup, 3)
        pileup = cleaner.run()

        correct_pileup = DensePileup.from_intervals(
            graph,
            [
                Interval(0, 6, [1])
            ]
        )

        self.assertEqual(pileup, correct_pileup)

    def test_single_hole_single_rp(self):
        pileup = DensePileup.from_intervals(graph,
                                            [Interval(0, 3, [1]),
                                             Interval(5, 7, [1])]
                                            )
        pileup.threshold(0.5)

        cleaner = DagHoleCleaner(pileup, 3)
        pileup = cleaner.run()

        correct_pileup = DensePileup.from_intervals(
            graph,
            [
                Interval(0, 10, [1])
            ]
        )

        self.assertEqual(pileup, correct_pileup)

    def test_single_hole_dual_rp(self):
        pileup = DensePileup.from_intervals(graph,
                                            [Interval(0, 8, [1]),
                                             Interval(3, 7, [2])]
                                            )
        pileup.threshold(0.5)

        cleaner = DagHoleCleaner(pileup, 5)
        pileup = cleaner.run()

        correct_pileup = DensePileup.from_intervals(
            graph,
            [
                Interval(0, 2, [1, 2, 3])
            ]
        )

        self.assertEqual(pileup, correct_pileup)

    def test_single_peak_split_graph(self):
        pileup = DensePileup.from_intervals(split_graph,
                                            [Interval(0, 10, [1])]
                                            )
        pileup.threshold(0.5)

        cleaner = DagHoleCleaner(pileup, 5)
        pileup = cleaner.run()

        correct_pileup = DensePileup.from_intervals(
            split_graph,
            [
                Interval(0, 5, [1, 2]),
                Interval(0, 5, [3])
            ]
        )

        self.assertEqual(pileup, correct_pileup)

    def test_single_hole_split_graph(self):
        pileup = DensePileup.from_intervals(split_graph,
                                            [Interval(0, 8, [1]),
                                             Interval(2, 5, [2])]
                                            )
        pileup.threshold(0.5)

        cleaner = DagHoleCleaner(pileup, 5)
        pileup = cleaner.run()

        correct_pileup = DensePileup.from_intervals(
            split_graph,
            [
                Interval(0, 10, [1, 2]),
                Interval(0, 3, [3])
            ]
        )
        print(pileup)

        self.assertEqual(pileup, correct_pileup)
"""


if __name__ == "__main__":
    unittest.main()