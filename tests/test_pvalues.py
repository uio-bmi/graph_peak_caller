from graph_peak_caller.pvalues import PToQValuesMapper, PValuesFinder, QValuesFinder
from graph_peak_caller.densepileup import DensePileup
from offsetbasedgraph import GraphWithReversals as Graph, DirectedInterval as Interval, Block
import unittest
import numpy as np

class TestPToQValuesMapper(unittest.TestCase):

    def test_simple_conversion(self):
        p_value_counts = {2.0: 2, 1: 1, 0.5: 3}
        mapper = PToQValuesMapper(p_value_counts)
        mapping = mapper.get_p_to_q_values()

        q_val_2 = 2 + (np.log10(1) - np.log10(6))
        self.assertAlmostEqual(mapping["%.7f" % 2], q_val_2)

        q_val_1 = 1 + (np.log10(3) - np.log10(6))
        self.assertAlmostEqual(mapping["%.7f" % 1], q_val_1)

        q_val_05 = 0.5 + (np.log10(4) - np.log10(6))
        self.assertAlmostEqual(mapping["%.7f" % 0.5], q_val_05)


class TestPValuesFinder(unittest.TestCase):


    def setUp(self):
        self.graph = Graph({i: Block(3) for i in range(1, 3)},
                  {1: [2]})

    def test_sample_equals_control(self):
        sample = DensePileup.from_intervals(self.graph, [Interval(0, 3, [1, 2])])
        control = DensePileup.from_intervals(self.graph, [Interval(0, 3, [1, 2])])

        finder = PValuesFinder(sample, control)
        p_values = finder.get_p_values_pileup()

        correct = -np.log10(0.26424)
        self.assertTrue(np.allclose(p_values.data._values, correct))

    def test_sample_equals_control_one_node(self):
        sample = DensePileup.from_intervals(self.graph, [Interval(0, 3, [2])])
        control = DensePileup.from_intervals(self.graph, [Interval(0, 3, [1, 2])])

        finder = PValuesFinder(sample, control)
        p_values = finder.get_p_values_pileup()

        correct = -np.log10(0.26424)
        self.assertTrue(np.allclose(p_values.data.values(1), 0))
        self.assertTrue(np.allclose(p_values.data.values(2), correct))

    def test_sample_twice_control_one_node(self):
        sample = DensePileup.from_intervals(self.graph, [Interval(0, 3, [2]), Interval(0, 3, [2])])
        control = DensePileup.from_intervals(self.graph, [Interval(0, 3, [1, 2])])

        finder = PValuesFinder(sample, control)
        p_values = finder.get_p_values_pileup()

        correct = -np.log10(0.08030)
        self.assertTrue(np.allclose(p_values.data.values(1), 0))
        self.assertTrue(np.allclose(p_values.data.values(2), correct))



if __name__ == "__main__":
    unittest.main()
