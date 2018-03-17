from graph_peak_caller.sparsepvalues import PToQValuesMapper, PValuesFinder
from graph_peak_caller.sparsediffs import SparseValues
from offsetbasedgraph import GraphWithReversals as Graph,\
    DirectedInterval as Interval, Block
import unittest
import numpy as np
from util import from_intervals


class TestPToQValuesMapper(unittest.TestCase):

    def test_simple_conversion(self):
        p_values = [2., 1., 0.5]
        counts = [2, 3, 6]
        mapper = PToQValuesMapper(p_values, counts)
        mapping = mapper.get_p_to_q_values()

        q_val_2 = 2 + (np.log10(1) - np.log10(6))
        self.assertAlmostEqual(mapping[2], q_val_2)

        q_val_1 = 1 + (np.log10(3) - np.log10(6))
        self.assertAlmostEqual(mapping[1], q_val_1)

        q_val_05 = 0.5 + (np.log10(4) - np.log10(6))
        self.assertAlmostEqual(mapping[0.5], q_val_05)


class TestPValuesFinder(unittest.TestCase):

    def setUp(self):
        self.graph = Graph({i: Block(3) for i in range(1, 3)},
                           {1: [2]})

    def test_sample_equals_control(self):
        sample = from_intervals(
            self.graph, [Interval(0, 3, [1, 2])])
        control = from_intervals(
            self.graph, [Interval(0, 3, [1, 2])])

        finder = PValuesFinder(sample, control)
        p_values = finder.get_p_values_pileup()
        correct = SparseValues([0, 6], [-np.log10(0.26424), 0])
        self.assertEqual(p_values, correct)
        # self.assertTrue(np.allclose(p_values.data._values, correct))

    def test_sample_equals_control_one_node(self):
        sample = from_intervals(
            self.graph, [Interval(0, 3, [2])])
        control = from_intervals(
            self.graph, [Interval(0, 3, [1, 2])])

        finder = PValuesFinder(sample, control)
        p_values = finder.get_p_values_pileup()
        correct = SparseValues([0, 3, 6], [0, -np.log10(0.26424), 0])
        self.assertEqual(p_values, correct)

    def test_sample_twice_control_one_node(self):
        sample = from_intervals(
            self.graph, [Interval(0, 3, [2]), Interval(0, 3, [2])])
        control = from_intervals(
            self.graph, [Interval(0, 3, [1, 2])])

        finder = PValuesFinder(sample, control)
        p_values = finder.get_p_values_pileup()
        correct = SparseValues([0, 3, 6], [0, -np.log10(0.08030), 0])
        self.assertEqual(p_values, correct)


if __name__ == "__main__":
    unittest.main()
