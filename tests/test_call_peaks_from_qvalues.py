import unittest
from offsetbasedgraph import Block, Interval, DirectedInterval, GraphWithReversals
from graph_peak_caller.callpeaks import CallPeaksFromQvalues, CallPeaks, ExperimentInfo
from graph_peak_caller.sparsepileup import SparsePileup, ValuedIndexes
import logging
logging.basicConfig(level=logging.ERROR)


class TestCallPeaksFromQValues(unittest.TestCase):

    def setUp(self):
        blocks = {i: Block(10) for i in range(1, 5)}
        edges = {i: [i+1] for i in range(1, 4)}
        self.linear_graph = GraphWithReversals(blocks, edges)
        self.one_peak_q_values = SparsePileup(self.linear_graph)
        self.one_peak_q_values.data = \
            {
                1: ValuedIndexes([5], [2], 0, 10),
                2: ValuedIndexes([3], [0], 2, 10)
            }

        self.one_peak_with_hole = SparsePileup(self.linear_graph)
        self.one_peak_with_hole.data = \
            {
                1: ValuedIndexes([5, 8], [2, 0], 0, 10),
                2: ValuedIndexes([3], [0], 2, 10)
            }

        self.fragment_length = 6
        self.read_length = 2


    def _run_caller(self, graph, pileup):
        graph_size = sum(block.length() for block in graph.blocks.values())
        experiment_info = ExperimentInfo(graph_size, self.fragment_length,
                                         self.read_length)
        caller = CallPeaksFromQvalues(graph, pileup, experiment_info,
                                    out_file_base_name="test_",
                                    cutoff=0.1)
        caller.callpeaks()
        return caller.max_paths

    def _assert_finds_max_paths(self, max_paths, graph, pileup):
        found_max_paths = self._run_caller(graph, pileup)
        for path in max_paths:
            self.assertTrue(path in found_max_paths)

        self.assertEqual(len(found_max_paths), len(max_paths))

    def test_finds_single_peak(self):
        self._assert_finds_max_paths([Interval(5, 3, [1, 2])],
                                     self.linear_graph, self.one_peak_q_values)

    def test_finds_single_peak_with_hole(self):
        self._assert_finds_max_paths([Interval(5, 3, [1, 2])],
                                     self.linear_graph, self.one_peak_with_hole)




if __name__ == "__main__":
    unittest.main()