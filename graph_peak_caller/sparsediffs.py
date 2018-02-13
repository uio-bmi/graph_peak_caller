import numpy as np


class SparseDiffs:
    def __init__(self, indices, diffs, sanitize=False):
        self._indices = np.asanyarray(indices)
        self._diffs = np.asanyarray(diffs)
        if sanitize:
            self._sanitize()

    def __repr__(self):
        return "SD(%s, %s)" % (self._indices, self._diffs)

    def maximum(self, other):
        all_indices = np.r_[self._indices, other._indices]
        sorted_args = np.argsort(all_indices, kind="mergesort")
        new_indexes = all_indices[sorted_args]
        unique_mask = np.ediff1d(new_indexes, to_end=1) != 0
        new_indexes = new_indexes[unique_mask]
        diffs = np.zeros((2, all_indices.size))
        my_args = sorted_args < self._indices.size
        other_args = ~my_args
        diffs[0][my_args] = self._diffs
        diffs[1][other_args] = other._diffs
        values = np.cumsum(diffs, axis=1)
        max_values = np.max(values, axis=0)
        max_values = max_values[unique_mask]
        new_diffs = np.ediff1d(
            max_values, to_begin=max_values[0])
        return SparseDiffs(new_indexes, new_diffs, True)

    def _sanitize(self):
        # Remove duplicate indexes
        index_diffs = np.ediff1d(self._indices, to_end=1)
        changes = index_diffs != 0
        self._indices = self._indices[changes]
        self._diffs = self._diffs[changes]

        # Remove duplicated values
        changes = self._diffs != 0
        self._indices = self._indices[changes]
        self._diffs = self._diffs[changes]

    def __eq__(self, other):
        return np.all(self._indices == other._indices)\
            and np.all(self._diffs == other._diffs)

    @classmethod
    def from_starts_and_ends(cls, starts_ends):
        indices = starts_ends.ravel(order="C")
        args = np.argsort(indices)
        indices = indices[args]
        diffs = np.ones(indices.size)
        diffs[args >= starts_ends.shape[1]] *= -1
        return cls(indices, diffs)

    def __itruediv__(self, scalar):
        self._diffs /= scalar
        return self


def _get_random_vec(n_points, size):
    from numpy.random import rand, randint
    return SparseDiffs(np.sort(randint(0, size, n_points)),
                       rand(n_points))

if __name__ == "__main__":
    import cProfile
    a = _get_random_vec(20000000, 500000000)
    b = _get_random_vec(20000000, 500000000)
    cProfile.run("a.maximum(b)")
