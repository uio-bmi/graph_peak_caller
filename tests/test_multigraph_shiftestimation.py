import unittest
from graph_peak_caller.shift_estimation_multigraph import MultiGraphShiftEstimator


class TestMultiGraphShiftEstimation(unittest.TestCase):

    def __test_simple(self):
        estimator = MultiGraphShiftEstimator(
            {"+":
                 {
                     1: [1, 11, 21],
                     2: [5, 15, 25]
                 },
             "-":
                 {
                    1:[6, 16, 26],
                    2: [10, 20, 30]
                }
            }
        )

        shift = estimator.get_estimates()


if __name__ == "__main__":
    unittest.main()








