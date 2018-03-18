# from graph_peak_caller.util import sparse_maximum
import unittest
import numpy as np


class _TestUtil:

    def test_sparse_maximum(self):

        indices1 = np.array([0, 4])
        values1 = np.array([0, 3])
        indices2 = np.array([0, 4])
        values2 = np.array([1, 2])

        max_indices, max_values = sparse_maximum(
            indices1, values1, indices2, values2, 10)
        self.assertTrue(np.all(max_indices == [0, 4]))
        self.assertTrue(np.all(max_values == [1, 3]))

    def test_sparse_maximum2(self):

        indices1 = np.array([1, 4])
        values1 = np.array([0, 3])
        indices2 = np.array([2, 4])
        values2 = np.array([1, 2])

        max_indices, max_values = sparse_maximum(indices1, values1,
                                                 indices2, values2, 10)
        print(max_indices)
        self.assertTrue(np.all(max_indices == [1, 2, 4]))
        self.assertTrue(np.all(max_values == [0, 1, 3]))

if __name__ == "__main__":
    unittest.main()
