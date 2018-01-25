import unittest
import numpy as np

from graph_peak_caller.linear_shift_estimation import LinearShiftEstimator




class TestLinearShiftEstimation(unittest.TestCase):

    def setUp(self):
        values = np.zeros(500)
        values[10:30] = 1
        values[100:150] += 1
        values[110:120] += 1
        values[200:300] += 1
        values[210:290] += 1
        values[250:260] += 1


        self.values = values

    def test_find_subpeaks(self):
        estimator = LinearShiftEstimator()
        subpeaks = estimator.find_subpeak_indices(self.values)

        self.assertTrue(np.all(subpeaks == [20, 115, 255]))

    def test_find_subpeak_pairs(self):
        estimator = LinearShiftEstimator()
        forward, reverse = estimator.find_subpeak_pairs(np.array([15]), np.array([25]))
        self.assertTrue(np.all(forward == [15]))
        self.assertTrue(np.all(reverse == [25]))

        forward, reverse = estimator.find_subpeak_pairs(np.array([15, 20]), np.array([21, 25]))
        self.assertTrue(np.all(forward == [20]))
        self.assertTrue(np.all(reverse == [21]))

    def test_full(self):
        forward = np.zeros(400)
        reverse = np.zeros(400)
        forward[10:20] += 1
        reverse[30:40] += 1

        forward[200:208] += 1
        reverse[210:218] += 1

        forward[300:308] += 1
        reverse[310:318] += 1

        estimator = LinearShiftEstimator(forward, reverse, min_pileup_value=0)
        fragment_length = estimator.run()
        self.assertEqual(fragment_length, 10)



if __name__ == "__main__":
    unittest.main()

