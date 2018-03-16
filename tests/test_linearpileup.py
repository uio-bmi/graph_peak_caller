import unittest
import numpy as np
from graph_peak_caller.control.linearsnarls import LinearPileup


class TestLinearPileup(unittest.TestCase):

    def setUp(self):
        starts = np.array([1, 4, 6])
        ends = np.array([2, 5, 10])
        self.pileup = LinearPileup.create_from_starts_and_ends(starts, ends)

        self.valued_pileup = LinearPileup(np.array([1, 4, 8]),
                                          np.array([1, 0, 3]))

    def test_create_from_starts_and_ends(self):
        self.assertEqual(LinearPileup(np.array([1, 2, 4, 5, 6, 10]),
                                      np.array([1., 0., 1., 0., 1., 0.])),
                         self.pileup)

    def test_maximum(self):
        pileup1 = LinearPileup(np.array([0,   3, 5.2, 10.8]),
                               np.array([10, 20,  15,   25]))

        pileup2 = LinearPileup(np.array([0,  5.2, 10]),
                               np.array([25,  10, 15]))
        pileup1.maximum(pileup2)

        true_max = LinearPileup(np.array([0, 5.2, 10.8]),
                                np.array([25, 15, 25]))
        self.assertEqual(pileup1, true_max)

if __name__ == "__main__":
    unittest.main()
