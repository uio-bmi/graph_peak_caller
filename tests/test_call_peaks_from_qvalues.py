import unittest
from offsetbasedgraph import Block, Interval, DirectedInterval, GraphWithReversals
from graph_peak_caller.callpeaks import CallPeaksFromQvalues, CallPeaks, ExperimentInfo
from graph_peak_caller.sparsepileup import SparsePileup, ValuedIndexes
from graph_peak_caller.densepileup import DensePileup, DensePileupData
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

        self.one_peak_with_big_hole = SparsePileup(self.linear_graph)
        self.one_peak_with_big_hole.data = \
            {
                1: ValuedIndexes([5, 7], [2, 0], 0, 10),
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

        self.single_block_graph = \
            GraphWithReversals({1: Block(20)}, {})

        self.multi_start_end_graph = \
                GraphWithReversals({i: Block(10) for i in range(1, 6)},
                           {
                               1: [3],
                               2: [3],
                               3: [4, 5]
                           })

        self.fragment_length = 6
        self.read_length = 2

    def _run_caller(self, graph, pileup):
        pileup = DensePileup.create_from_old_sparsepileup(pileup)
        graph_size = sum(block.length() for block in graph.blocks.values())
        experiment_info = ExperimentInfo(graph_size, self.fragment_length,
                                         self.read_length)
        caller = CallPeaksFromQvalues(graph, pileup, experiment_info,
                                      out_file_base_name="test_",
                                      cutoff=0.1, q_values_max_path=True)
        caller.callpeaks()
        return caller.max_paths

    def _assert_finds_max_paths(self, max_paths, graph, pileup):
        found_max_paths = self._run_caller(graph, pileup)
        for path in max_paths:
            path.graph = graph
            self.assertTrue(path in found_max_paths or path.get_reverse() in found_max_paths,
                "\nPath %s not found in max paths. Max paths: \n %s" % \
                (path, '\n'.join([str(p) for p in found_max_paths])))

        self.assertEqual(len(found_max_paths), len(max_paths))

    def test_finds_single_peak(self):
        self._assert_finds_max_paths([Interval(5, 3, [1, 2])],
                                     self.linear_graph, self.one_peak_q_values)

    def test_finds_single_peak_with_hole(self):
        self._assert_finds_max_paths([Interval(5, 3, [1, 2])],
                                     self.linear_graph, self.one_peak_with_hole)

    def test_finds_no_peak_too_large_hole(self):
        self._assert_finds_max_paths([],
                                     self.linear_graph, self.one_peak_with_big_hole)

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

    def test_removes_too_small_peak(self):
        pileup = SparsePileup(self.single_block_graph)
        pileup.data = {
            1: ValuedIndexes([5, 9], [2, 0], 0, 20)
        }
        self._assert_finds_max_paths([],
                    self.single_block_graph,
                    pileup)

    def test_finds_internal_peak_with_internal_hole(self):
        pileup = SparsePileup(self.single_block_graph)
        pileup.data = {
            1: ValuedIndexes([5, 8, 9, 15], [2, 0, 2, 0], 0, 20)
        }
        self._assert_finds_max_paths(
                    [Interval(5, 15, [1])],
                    self.single_block_graph,
                    pileup)

    def test_finds_correct_max_path_among_many_paths(self):
        graph = GraphWithReversals(
            {
                1: Block(10),
                2: Block(10),
                3: Block(10),
                4: Block(10),
                5: Block(10)
            },
            {
                1: [2, 3, 4],
                2: [5],
                4: [5],
                4: [5]
            }
        )

        pileup = SparsePileup(graph)
        pileup.data = {
            1: ValuedIndexes([], [], 2, 10),
            # Higher qval, but two holes with low
            2: ValuedIndexes([1, 2, 7, 8], [0, 2.001, 0, 2.001], 2, 10),
            3: ValuedIndexes([], [], 1.5, 10),
            4: ValuedIndexes([], [], 2, 10),
            5: ValuedIndexes([], [], 2, 10)
        }
        self._assert_finds_max_paths(
            [Interval(0, 10, [1, 4, 5])],
            graph, pileup
        )

    def test_multiple_start_and_end_nodes(self):

        pileup = SparsePileup(self.multi_start_end_graph)
        pileup.data = {
            1: ValuedIndexes([], [], 2, 10),
            2: ValuedIndexes([], [], 2.2, 10),
            3: ValuedIndexes([1, 9], [2, 0], 0, 10),
            4: ValuedIndexes([], [], 2, 10),
            5: ValuedIndexes([3], [3], 0, 10),
        }

        self._assert_finds_max_paths(
            [
                Interval(0, 10, [2, 3, 4]),
                Interval(3, 10, [5])
            ],
            self.multi_start_end_graph, pileup
        )

    def test_multiple_start_and_end_nodes2(self):
        pileup = SparsePileup(self.multi_start_end_graph)
        pileup.data = {
            1: ValuedIndexes([], [], 2, 10),
            2: ValuedIndexes([], [], 2.2, 10),
            3: ValuedIndexes([1, 9], [2, 0], 0, 10),
            4: ValuedIndexes([], [], 2, 10),
            5: ValuedIndexes([1], [3], 0, 10),
        }

        self._assert_finds_max_paths(
            [
                Interval(0, 10, [2, 3, 5])
            ],
            self.multi_start_end_graph, pileup
        )


if __name__ == "__main__":
    unittest.main()
