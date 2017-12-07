import unittest
from offsetbasedgraph import GraphWithReversals, Block, \
        Interval, DirectedInterval, IntervalCollection
from graph_peak_caller.callpeaks import ExperimentInfo, CallPeaks
from graph_peak_caller.snarls import SnarlGraph, SnarlGraphBuilder, SimpleSnarl
from graph_peak_caller.linearsnarls import LinearSnarlMap, create_control
from graph_peak_caller.sparsepileup import SparsePileup, ValuedIndexes
import numpy as np


class TestCreateControl(unittest.TestCase):

    def setUp(self):
        self.graph = GraphWithReversals(
            {i: Block(3) for i in range(1, 13)},
            {
                11: [1],
                1: [2, 3],
                2: [7, 8],
                3: [4, 5],
                4: [6],
                5: [6],
                6: [10],
                7: [9],
                8: [9],
                9: [10],
                10: [12]
             })

        self.linear_length = 21
        self.snarlgraph = SnarlGraph(
            {
                11: Block(3),
                12: Block(3),
                1: Block(3),
                10: Block(3),
                20: SnarlGraph(
                    {
                        3: Block(3),
                        21: SnarlGraph(
                            {
                                4: Block(3),
                                5: Block(3)
                            },
                            {
                                3: [4, 5],
                                4: [6],
                                5: [6]
                            },
                            start_node=3,
                            end_node=6
                        ),
                        22: SnarlGraph(
                            {
                                7: Block(3),
                                8: Block(3)
                            },
                            {
                                2: [7, 8],
                                7: [9],
                                8: [9]
                            },
                            start_node=2,
                            end_node=9
                        ),
                        2: Block(3),
                        6: Block(3),
                        9: Block(3),
                    },
                    {
                        3: [21],
                        2: [22],
                        21: [6],
                        22: [9],
                        1: [2, 3],
                        6: [10],
                        9: [10]
                    },
                    start_node=1,
                    end_node=10
                )
            },
            {
                11: [1],
                1: [20],
                20: [10],
                10: [12],
                13: [11],  # Dummy
                12: [14],   # Dummy
            },
            start_node=13,
            end_node=14
        )

        LinearSnarlMap.from_snarl_graph(
            self.snarlgraph, self.graph).to_json_files("test_linear_map.tmp")

    def test_single_read_single_extension(self):
        fragment_length = 3
        reads = [Interval(0, 3, [2])]
        extension_sizes = [8]
        control = create_control("test_linear_map.tmp", reads, extension_sizes,
                                 fragment_length, ob_graph=self.graph)

        expected_bakground = len(reads) * fragment_length / self.linear_length
        value_in_extension = 1 * fragment_length / (extension_sizes[0])

        correct_pileup = SparsePileup.from_base_value(self.graph, expected_bakground)
        for rp in [2, 3, 1]:
            correct_pileup.data[rp] = ValuedIndexes([], [], value_in_extension, 3)

        for rp in [7, 8, 4, 5]:
            correct_pileup.data[rp] = ValuedIndexes([1], [expected_bakground], value_in_extension, 3)

        for rp in [11]:
            correct_pileup.data[rp] = ValuedIndexes([2], [value_in_extension], expected_bakground, 3)

        self.assertEqual(control, correct_pileup)

    def test_single_read_two_extensions(self):
        fragment_length = 3
        reads = [Interval(0, 3, [2])]
        extension_sizes = [2, 8]
        control = create_control("test_linear_map.tmp", reads, extension_sizes, fragment_length, ob_graph=self.graph)

        expected_bakground = len(reads) * fragment_length / self.linear_length
        value_in_extensions = 1 * fragment_length / (np.array(extension_sizes))

        correct_pileup = SparsePileup.from_base_value(self.graph, expected_bakground)
        for rp in [2, 3]:
            correct_pileup.data[rp] = ValuedIndexes([1], [value_in_extensions[1]], value_in_extensions[0], 3)

        for rp in [1]:
            correct_pileup.data[rp] = ValuedIndexes([2], [value_in_extensions[0]], value_in_extensions[1], 3)

        for rp in [7, 8, 4, 5]:
            correct_pileup.data[rp] = ValuedIndexes([1], [expected_bakground], value_in_extensions[1], 3)

        for rp in [11]:
            correct_pileup.data[rp] = ValuedIndexes([2], [value_in_extensions[1]], expected_bakground, 3)

        self.assertEqual(control, correct_pileup)

    def test_two_reads_two_extensions(self):
        fragment_length = 3
        reads = [Interval(0, 3, [2]), Interval(1, 1, [8, 9])]
        extension_sizes = [2, 8]
        control = create_control("test_linear_map.tmp", reads, extension_sizes, fragment_length, ob_graph=self.graph)

        expected_bakground = len(reads) * fragment_length / self.linear_length
        value_in_extensions = 1 * fragment_length / (np.array(extension_sizes))
        correct_pileup = SparsePileup.from_base_value(self.graph, expected_bakground)
        for rp in [2, 3]:
            correct_pileup.data[rp] = ValuedIndexes([1], [value_in_extensions[1]*2], value_in_extensions[0], 3)

        for rp in [1]:
            correct_pileup.data[rp] = ValuedIndexes([2], [value_in_extensions[0]], value_in_extensions[1], 3)

        for rp in [7, 8, 4, 5]:
            correct_pileup.data[rp] = ValuedIndexes([2], [value_in_extensions[1]], value_in_extensions[0], 3)

        for rp in [6, 9]:
            correct_pileup.data[rp] = ValuedIndexes([2], [expected_bakground], value_in_extensions[1], 3)

        for rp in [11]:
            correct_pileup.data[rp] = ValuedIndexes([2], [value_in_extensions[1]], expected_bakground, 3)

        #print(correct_pileup)
        #print(control)
        for i in range(1, 13):
            self.assertEqual(control.data[i], correct_pileup.data[i], "%d not equal" % i)

        self.assertEqual(control, correct_pileup)


class _TestCreateControlGraphWithDifferentLengths():

    def setUp(self):
        self.graph = GraphWithReversals(
            {i: Block(3) for i in range(1, 13)},
            {
                11: [1],
                1: [2, 3],
                2: [7, 8],
                3: [6],
                6: [10],
                7: [9],
                8: [9],
                9: [10],
                10: [12]
             })
        self.blocks[3] = Block(6)

        self.linear_length = 21
        self.snarlgraph = SnarlGraph(
            {
                11: Block(3),
                12: Block(3),
                1: Block(3),
                10: Block(3),
                20: SnarlGraph(
                    {
                        3: Block(6),
                        22: SnarlGraph(
                            {
                                7: Block(3),
                                8: Block(3)
                            },
                            {
                                2: [7, 8],
                                7: [9],
                                8: [9]
                            },
                            start_node=2,
                            end_node=9
                        ),
                        2: Block(3),
                        6: Block(3),
                        9: Block(3),
                    },
                    {
                        3: [6],
                        2: [22],
                        22: [9],
                        1: [2, 3],
                        6: [10],
                        9: [10]
                    },
                    start_node=1,
                    end_node=10
                )
            },
            {
                11: [1],
                1: [20],
                20: [10],
                10: [12],
                13: [11],  # Dummy
                12: [14],   # Dummy
            },
            start_node=13,
            end_node=14
        )

        LinearSnarlMap.from_snarl_graph(self.snarlgraph, self.graph).to_json_files("test_linear_map.tmp")

    def test_single_read(self):
        fragment_length = 3
        reads = [Interval(0, 3, [2])]
        extension_sizes = [8]
        control = create_control("test_linear_map.tmp", reads, extension_sizes, fragment_length, ob_graph=self.graph)
        expected_bakground = len(reads) * fragment_length / self.linear_length
        value_in_extension = 1 * fragment_length / (extension_sizes[0])

        correct_pileup = SparsePileup.from_base_value(self.graph, expected_bakground)
        for rp in [2, 3, 1]:
            correct_pileup.data[rp] = ValuedIndexes([], [], value_in_extension, 3)

        for rp in [7, 8, 4, 5]:
            correct_pileup.data[rp] = ValuedIndexes([1], [expected_bakground], value_in_extension, 3)

        for rp in [11]:
            correct_pileup.data[rp] = ValuedIndexes([2], [value_in_extension], expected_bakground, 3)

        self.assertEqual(control, correct_pileup)

if __name__ == "__main__":
    unittest.main()
