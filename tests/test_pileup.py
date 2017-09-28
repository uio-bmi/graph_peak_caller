import unittest
from graph_peak_caller.pileup import Pileup, SmallIntervals

from examples import *


class TestPileup(unittest.TestCase):
    def test_create(self):
        pileup = Pileup(ob_graphs[0])
        pileup.add_intervals(directed_pileup_intervals)
        self.assertEqual(pileup, true_pileup)

    def test_to_file_from_file(self):
        pileup =  Pileup(ob_graphs[0])
        pileup.add_intervals(pileup_intervals)
        pileup.to_bed_graph("test.pileup")

        pileup_from_file = Pileup.from_bed_graph(ob_graphs[0], "test.pileup")

        self.assertTrue(pileup_from_file == pileup)

    def test_fill_small_holes(self):
        graph = ob_graphs[0]
        init_counts = {node.id: np.zeros(node.n_basepairs, dtype="int32")
                       for node in nodes}
        init_counts[0][3:8] = 1
        init_counts[3][3:8] = 1
        pileup = Pileup(graph)
        pileup.set_count_arrays(init_counts)
        pileup.fill_small_wholes(26)
        true_counts = {node.id: np.zeros(node.n_basepairs, dtype="int32")
                       for node in nodes}
        true_counts[0][3:10] = 1
        true_counts[3][0:8] = 1
        true_counts[1][0:20] = 1
        true_pileup = Pileup(graph)
        true_pileup.set_count_arrays(true_counts)
        self.assertEqual(pileup, true_pileup)


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
