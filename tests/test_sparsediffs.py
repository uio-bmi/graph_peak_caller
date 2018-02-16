import unittest
import numpy as np
from graph_peak_caller.sparsediffs import SparseDiffs, SparseValues


class TestSparseValues(unittest.TestCase):

    def test_to_from_file(self):

        indexes = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        values = np.array([1, 2, 3, 4, 5, 6, 7, 8])

        sv = SparseValues(indexes, values)
        sv.track_size = 10
        sv.to_sparse_files("test_sparsevalues.tmp")

        new = sv.from_sparse_files("test_sparsevalues.tmp")
        self.assertEqual(sv, new)


if __name__ == "__main__":
    unittest.main()