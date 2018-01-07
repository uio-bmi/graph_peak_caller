from graph_peak_caller.sparsepileupv2 import SparsePileup, SparsePileupData, SparseControlSample
from graph_peak_caller.sparsepileup import SparsePileup as SparsePileupOld
from offsetbasedgraph import GraphWithReversals, Block, Interval

import unittest
import numpy as np

class TestSparsePileup(unittest.TestCase):

    def test_from_starts_and_ends(self):

        graph = GraphWithReversals({i: Block(10) for i in range(1, 4)},
                                   {1: [2],
                                    2: [3]})

        starts = {1: [3, 5]}
        ends = {1: [7, 9]}
        pileup = SparsePileup.from_starts_and_ends(graph, starts, ends)
        print(pileup)
        self.assertTrue(np.all(pileup.data.values(1) == [0, 1, 2, 1, 0]))
        self.assertTrue(np.all(pileup.data.indexes(1) == [0, 3, 5, 7, 9, 10]))

        starts = {1: [0, 3]}
        ends = {1: [5, 10]}
        pileup = SparsePileup.from_starts_and_ends(graph, starts, ends)
        self.assertTrue(np.all(pileup.data.values(1) == [1, 2, 1]))
        self.assertTrue(np.all(pileup.data.indexes(1) == [0, 3, 5, 10]))

        starts = {1: [0, 3]}
        ends = {1: [5, 9]}
        pileup = SparsePileup.from_starts_and_ends(graph, starts, ends)
        self.assertTrue(np.all(pileup.data.values(1) == [1, 2, 1, 0]))
        self.assertTrue(np.all(pileup.data.indexes(1) == [0, 3, 5, 9, 10]))

        starts = {1: [0], 2: [0]}
        ends = {1: [10], 2: [10]}
        pileup = SparsePileup.from_starts_and_ends(graph, starts, ends)

        self.assertTrue(np.all(pileup.data.values(1) == [1]))
        self.assertTrue(np.all(pileup.data.indexes(1) == [0, 10]))

        self.assertTrue(np.all(pileup.data.values(2) == [1]))
        self.assertTrue(np.all(pileup.data.indexes(2) == [0, 10]))


        """
        old = SparsePileupOld.from_starts_and_ends(graph, starts, ends)
        print(old)
        print(old.data[1].all_values())
        print(old.data[1].all_idxs())
        """


class TestSparseControlSample(unittest.TestCase):

    def test_from_control_and_sample(self):
        graph = GraphWithReversals({i: Block(10) for i in range(1, 4)},
                                   {1: [2],
                                    2: [3]})

        control = SparsePileup.from_starts_and_ends(
            graph,
            {1: [0, 4]},
            {1: [5, 10]},
        )
        sample = SparsePileup.from_starts_and_ends(
            graph,
            {1: [0, 7]},
            {1: [8, 10]},
        )

        joined = SparseControlSample.from_sparse_control_and_sample(control, sample)
        self.assertTrue(np.all([0, 4, 5, 7, 8, 10] == joined.data.indexes(1)))
        self.assertTrue(np.all([1, 2, 1, 1, 1] == joined.data.values(1)[:,0]))
        self.assertTrue(np.all([1, 1, 1, 2, 1] == joined.data.values(1)[:,1]))

        control = SparsePileup.from_starts_and_ends(
            graph,
            {1: [0, 4]},
            {1: [5, 9]},
        )
        sample = SparsePileup.from_starts_and_ends(
            graph,
            {1: [0, 7]},
            {1: [8, 9]},
        )

        joined = SparseControlSample.from_sparse_control_and_sample(control, sample)
        self.assertTrue(np.all([0, 4, 5, 7, 8, 9, 10] == joined.data.indexes(1)))
        self.assertTrue(np.all([1, 2, 1, 1, 1, 0] == joined.data.values(1)[:,0]))
        self.assertTrue(np.all([1, 1, 1, 2, 1, 0] == joined.data.values(1)[:,1]))


if __name__ == "__main__":
    unittest.main()