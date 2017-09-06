import unittest
import numpy as np
from graph_peak_caller.sparsepileup import ValuedIndexes, \
    SparsePileup, SparseControlSample
from examples import *
from scipy.stats import poisson


def valued_indexes():
    return ValuedIndexes(np.array([10, 15], dtype="int"),
                         np.array([100, 200]),
                         50,
                         20)


class TestValuedIndexes(unittest.TestCase):
    def test_threshold(self):
        true_indexes = ValuedIndexes(np.array([15], dtype="int"),
                                     np.array([True], dtype="bool"),
                                     False,
                                     20)
        vi = valued_indexes()
        vi.threshold(150)
        vi.sanitize()
        self.assertEqual(vi, true_indexes)

    def test_set_interval_value_mid(self):
        vi = valued_indexes()
        vi.set_interval_value(10, 15, 150)
        true_vi = valued_indexes()
        true_vi.values[0] = 150
        self.assertEqual(vi, true_vi)

    def test_set_interval_value_start(self):
        vi = valued_indexes()
        vi.set_interval_value(0, 10, 150)
        true_vi = valued_indexes()
        true_vi.start_value = 150
        self.assertEqual(vi, true_vi)

    def test_set_interval_value_end(self):
        vi = valued_indexes()
        vi.set_interval_value(15, 20, 150)
        true_vi = valued_indexes()
        true_vi.values[1] = 150
        self.assertEqual(vi, true_vi)


class TestSparsePileup(unittest.TestCase):
    pass


class TestSparseControlSample(unittest.TestCase):
    def setUp(self):
        control = pileup1_one_block
        sample = pileup2_one_block
        a = control.get_count_array(1)
        a[0] = 3
        self.sparse = SparseControlSample.from_control_and_sample(
            control, sample)

    def test_from_control_and_sample(self):
        true_vi = ValuedIndexes(np.array([1, 5]),
                                np.array([[2., 1.],
                                          [1., 0.]]),
                                np.array([3., 1.]),
                                10)
        self.assertEqual(self.sparse.data[1], true_vi)

    def test_get_p_dict(self):
        self.sparse.get_p_dict()
        keys = [(2., 1.), (1., 0.), (3., 1.)]
        s_keys = [str(np.array(k)) for k in keys]
        p_values = [-np.log10(1-poisson.cdf(k[1], k[0])) for k in keys]
        true_dict = dict(zip(s_keys, p_values))
        self.assertEqual(self.sparse.p_value_dict,
                         true_dict)

        lens = [4, 5, 1]
        count_dict = dict(zip(p_values, lens))
        self.assertEqual(self.sparse.count_dict,
                         count_dict)

        self.sparse.get_p_to_q_values()
        print("#", self.sparse.p_to_q_values)
        self.sparse.get_q_values()

if __name__ == "__main__":
    unittest.main()
