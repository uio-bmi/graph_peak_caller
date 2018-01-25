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

    def test_to_from_sparse_files(self):
        pileup = DensePileup.from_intervals(graph,
                                            [
                                                Interval(5, 5, [1, 2], graph),
                                                Interval(7, 3, [1, 2], graph)
                                            ])

        pileup.to_sparse_files("test_pileup_")
        new = pileup.from_sparse_files(graph, "test_pileup_")

        self.assertEqual(new, pileup)
        self.assertEqual(pileup, new)

        pileup = DensePileup.from_intervals(graph,
                                            [
                                                Interval(0, 5, [1, 2], graph),
                                                Interval(7, 10, [1, 2], graph)
                                            ])

        pileup.to_sparse_files("test_pileup_")
        new = pileup.from_sparse_files(graph, "test_pileup_")

        self.assertEqual(new, pileup)
        self.assertEqual(pileup, new)

    def test_to_from_sparse_files2(self):
        pileup = DensePileup(graph)
        pileup.data.set_values(1, 1, 3, [3.5, 2.1])
        pileup.data.set_values(2, 0, 3, [1.1, 3.5, 2.1])
        pileup.data.set_values(3, 0, 4, [3.5, 3.5, 2.1, 2.1])

        pileup.to_sparse_files("test_pileup2_")
        new = DensePileup.from_sparse_files(graph, "test_pileup2_")
        self.assertTrue(np.all(new.data.values_in_range(1, 1, 3) == [3.5, 2.1]))
        self.assertTrue(np.all(new.data.values_in_range(2, 0, 3) == [1.1, 3.5, 2.1]))
        self.assertTrue(np.all(new.data.values_in_range(3, 0, 4) == [3.5, 3.5, 2.1, 2.1]))

    def test_to_bedgraph(self):
        pileup = DensePileup(graph)
        pileup.data.set_values(1, 3, 6, [1, 1, 1])
        pileup.data.set_values(1, 6, 8, [1.5, 2])

        pileup.to_bed_graph("test_bdg.tmp")

        lines = []
        f = open("test_bdg.tmp")
        for line in f:
            lines.append(line)
        f.close()

        self.assertTrue("1\t0\t3\t0.0\n" in lines)
        self.assertTrue("1\t3\t6\t1.0\n" in lines)
        self.assertTrue("1\t6\t7\t1.5\n" in lines)
        self.assertTrue("1\t7\t8\t2.0\n" in lines)

    def _test_to_from_bedgraph_equals_old_sparse_pileup(self):
        # Not working.
        pileup = DensePileup(graph)
        pileup.data.set_values(1, 3, 6, [1, 1, 1])
        pileup.data.set_values(2, 6, 8, [1.5, 2])
        pileup.to_bed_graph("test_bdg.tmp")
        from graph_peak_caller.sparsepileup import SparsePileup
        sparse = SparsePileup.from_bed_graph(graph, "test_bdg.tmp")
        self.assertTrue(pileup.equals_old_sparse_pileup(sparse))


if __name__ == "__main__":
    unittest.main()