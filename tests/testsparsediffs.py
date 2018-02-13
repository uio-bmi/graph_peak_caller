from graph_peak_caller.sparsediffs import SparseDiffs
import numpy as np


def test_maximum():
    a = SparseDiffs([0, 1, 3, 5, 7, 10], [1, 1, 1, -1, -1, -1])
    b = SparseDiffs([0, 2, 3, 6, 9, 10], [2, -1, 1, 1, -1, -2])
    assert a.maximum(b) == SparseDiffs([0, 3, 5, 6, 9, 10],
                                       [2, 1, -1, 1, -1, -2])


def test_from_startends():
    a = np.array([[1, 10, 9, 7, 11],
                  [3, 12, 11, 9, 13]])
    vec = SparseDiffs.from_starts_and_ends(a)
    assert vec == SparseDiffs([1, 3, 7, 9, 9, 10, 11, 11, 12, 13],
                              [1, -1, 1, 1, -1, 1, 1, -1, -1, -1])
