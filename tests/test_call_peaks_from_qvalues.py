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


        self.split_graph = GraphWithReversals(
            {i: Block(10) for i in range(1, 5)},
            {
                1: [2, 3],
                2: [4],
                3: [4]
            }
        )

        self.split_graph_with_path_around = \
            GraphWithReversals(
                {i: Block(10) for i in range(1, 8)},
                {}
            )

        self.graph_with_reversal = \
            GraphWithReversals(
                {i: Block(10) for i in range(1, 4)},
                {
                    1: [2],
                    -3: [-2]
                }
            )

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
            self.assertTrue(path in found_max_paths,
                "\nPath %s not found in max paths. Max paths: \n %s" % \
                (path, '\n'.join([str(p) for p in found_max_paths])))

        self.assertEqual(len(found_max_paths), len(max_paths))

    def test_finds_single_peak(self):
        self._assert_finds_max_paths([Interval(5, 3, [1, 2])],
                                     self.linear_graph, self.one_peak_q_values)

    def test_finds_single_peak_with_hole(self):
        self._assert_finds_max_paths([Interval(5, 3, [1, 2])],
                                     self.linear_graph, self.one_peak_with_hole)

    def test_find_max_path_on_split_graph(self):

        pileup = SparsePileup(self.split_graph)
        pileup.data = {
            1: ValuedIndexes([], [], 2, 10),
            2: ValuedIndexes([], [], 3, 10),
            3: ValuedIndexes([], [], 2, 10),
            4: ValuedIndexes([1, 4], [0, 3], 2, 10)
        }
        self._assert_finds_max_paths(
            [Interval(0, 1, [1, 2, 4]),
             Interval(4, 10, [4])],
            self.split_graph, pileup
        )

    def test_find_max_path_on_split_graph_with_hole_fill(self):

        pileup = SparsePileup(self.split_graph)
        pileup.data = {
            1: ValuedIndexes([3], [2], 0, 10),
            2: ValuedIndexes([], [], 2.1, 10),
            3: ValuedIndexes([], [], 2, 10),
            4: ValuedIndexes([1, 3], [0, 3], 2, 10)
        }
        self._assert_finds_max_paths(
            [Interval(3, 10, [1, 2, 4])],
            self.split_graph, pileup
        )

    def test_with_reversal(self):
        pileup = SparsePileup(self.graph_with_reversal)
        pileup.data = {
            1: ValuedIndexes([], [], 2, 10),
            2: ValuedIndexes([], [], 2, 10),
            3: ValuedIndexes([], [], 2.1, 10),
        }
        self._assert_finds_max_paths(
            [Interval(0, 10, [-3, -2])],
            self.graph_with_reversal,
            pileup
        )

    def test_with_reversal_and_hole(self):
        pileup = SparsePileup(self.graph_with_reversal)
        pileup.data = {
            1: ValuedIndexes([], [], 2, 10),
            2: ValuedIndexes([9], [0], 2, 10),
            3: ValuedIndexes([1], [3], 0, 10),
        }
        self._assert_finds_max_paths(
            [Interval(0, 10, [-3, -2])],
            self.graph_with_reversal,
            pileup
        )





if __name__ == "__main__":
    unittest.main()