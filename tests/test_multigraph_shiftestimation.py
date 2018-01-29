import unittest
from graph_peak_caller.shift_estimation_multigraph import MultiGraphShiftEstimator
import numpy as np
import logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s, %(levelname)s: %(message)s")

class TestMultiGraphShiftEstimation(unittest.TestCase):

    def test_simple(self):

        shift = 150
        genome_size = 10000000
        n_peaks = 5000

        positions = {"+":
                   {1: []},
               "-":
                   {1: []}
               }

        for i in range(0, n_peaks):
            pos = np.random.randint(0, genome_size-shift)

            for j in range(0, 40):
                pos += np.random.randint(-3, 3)
                positions["+"][1].append(pos)
                positions["-"][1].append(pos + shift)

        estimator = MultiGraphShiftEstimator(
            positions, genome_size
        )

        found_shift = estimator.get_estimates()
        self.assertEqual(round(found_shift), shift)
        print(shift)

if __name__ == "__main__":
    unittest.main()








