from graph_peak_caller.densepileup import DensePileupData, DensePileup
from offsetbasedgraph import GraphWithReversals as Graph, Block, DirectedInterval as Interval
import unittest
import numpy as np

graph = Graph({i: Block(10) for i in range(1, 4)},
                                   {1: [2],
                                    2: [3]})

class TestDensePileup(unittest.TestCase):


    def test_init(self):
        pileup = DensePileup(graph)

    def test_from_starts_and_ends(self):

        starts = {1: [3, 5]}
        ends = {1: [7, 9]}
        pileup = DensePileup.from_starts_and_ends(graph, starts, ends)
        indexes, values = pileup.data.get_sparse_indexes_and_values(1)
        self.assertTrue(np.all(values == [0, 1, 2, 1, 0]))
        self.assertTrue(np.all(indexes == [0, 3, 5, 7, 9, 10]))

        starts = {1: [0, 3]}
        ends = {1: [5, 10]}
        pileup = DensePileup.from_starts_and_ends(graph, starts, ends)
        indexes, values = pileup.data.get_sparse_indexes_and_values(1)
        self.assertTrue(np.all(values == [1, 2, 1]))
        self.assertTrue(np.all(indexes == [0, 3, 5, 10]))

        starts = {1: [0, 3]}
        ends = {1: [5, 9]}
        pileup = DensePileup.from_starts_and_ends(graph, starts, ends)
        indexes, values = pileup.data.get_sparse_indexes_and_values(1)
        self.assertTrue(np.all(values == [1, 2, 1, 0]))
        self.assertTrue(np.all(indexes == [0, 3, 5, 9, 10]))

        starts = {1: [0], 2: [0]}
        ends = {1: [10], 2: [10]}
        pileup = DensePileup.from_starts_and_ends(graph, starts, ends)
        indexes, values = pileup.data.get_sparse_indexes_and_values(1)

        self.assertTrue(np.all(values == [1]))
        self.assertTrue(np.all(indexes == [0, 10]))

        indexes, values = pileup.data.get_sparse_indexes_and_values(2)
        self.assertTrue(np.all(values == [1]))
        self.assertTrue(np.all(indexes == [0, 10]))


    def test_indexes_to_positions(self):
        starts = {1: [3, 5]}
        ends = {1: [7, 9]}
        pileup = DensePileup.from_starts_and_ends(graph, starts, ends)

        positions = pileup.data.value_indexes_to_nodes_and_offsets([0, 1, 10, 11, 15])
        self.assertEqual(positions, [(1, 0), (1, 1), (2, 0), (2, 1), (2, 5)])

        positions = pileup.data.value_indexes_to_nodes_and_offsets([15])
        self.assertEqual(positions, [(2, 5)])

        positions = pileup.data.value_indexes_to_nodes_and_offsets([15, 25])
        self.assertEqual(positions, [(2, 5), (3, 5)])

    def test_get_interval_values(self):
        pileup = DensePileup.from_intervals(graph,
                                            [
                                                Interval(5, 5, [1, 2], graph),
                                                Interval(7, 3, [1, 2], graph)
                                            ])
        values = pileup.data.get_interval_values(Interval(5, 5, [1, 2], graph))
        self.assertTrue(np.all(values == [1, 1, 2, 2, 2, 2, 2, 2, 1, 1]))

        values = pileup.data.get_interval_values(Interval(3, 6, [1], graph))
        self.assertTrue(np.all(values == [0, 0, 1]))


if __name__ == "__main__":
    unittest.main()