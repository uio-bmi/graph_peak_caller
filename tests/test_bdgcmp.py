import unittest
from graph_peak_caller.pileup import Pileup
from offsetbasedgraph import Graph, Block, Interval
import offsetbasedgraph
from graph_peak_caller import bdgcmp
import numpy as np
from examples import *
from graph_peak_caller.bdgcmp import *



class TestPileup(unittest.TestCase):
    def _test_simple_bdgcmp(self):
        graph = one_block_graph
        pileup1 = pileup1_one_block
        pileup2 = pileup2_one_block

        max_pileup = bdgcmp.create_background_pileup_as_max_from_pileups(graph, iter([pileup1, pileup2]), 0)

        correct_count_array = np.array([1, 2, 2, 2, 2, 1, 1, 1, 1, 1])
        closeness = np.isclose(max_pileup.get_count_arrays()[1], correct_count_array)
        print(closeness)
        print(max_pileup.get_count_arrays()[1])
        self.assertTrue(np.all(np.isclose(max_pileup.get_count_arrays()[1], correct_count_array)))

    def _test_get_p_value_track(self, graph, sample, control):

        pval_pileup_macs = get_p_value_track(graph, control.to_bed_graph("pileup2.tmp"), sample.to_bed_graph("pileup1.tmp"), "out.tmp")
        print("pval_pileup")
        print(pval_pileup_macs)

        pval_pileup_us = get_p_value_track_from_pileups(graph, control, sample)
        print("Pval pileup us")
        print(pval_pileup_us)
        pval_pileup_us.to_bed_graph("us.tmp")

        self.assertTrue(pval_pileup_macs.is_numerically_equal(pval_pileup_us))

    def test_get_p_value_track(self):
        graph = one_block_graph

        sample_intervals = [Interval(1, 10, [1], graph),
                            Interval(5, 7, [1], graph)]
        control_intervals = [Interval(0, 10, [1], graph)]
        sample = Pileup(graph)
        sample.add_intervals(sample_intervals)
        control = Pileup(graph)
        control.add_intervals(control_intervals)

        self._test_get_p_value_track(graph, sample, control)

        # case 2
        graph = offsetbasedgraph.Graph({1: Block(10), 2: Block(10)}, {1: [2]})
        sample_intervals = [Interval(1, 10, [1], graph),
                            Interval(5, 7, [1], graph),
                            Interval(4, 6, [1], graph),
                            Interval(0, 10, [2], graph),
                            Interval(1, 3, [2], graph)
                            ]

        control_intervals = [
                                Interval(0, 10, [1], graph),
                                Interval(0, 10, [2], graph),
                                Interval(5, 5, [1, 2], graph)
                             ]

        sample = Pileup(graph)
        sample.add_intervals(sample_intervals)
        control = Pileup(graph)
        control.add_intervals(control_intervals)

        self._test_get_p_value_track(graph, sample, control)

if __name__ == "__main__":
    unittest.main()
