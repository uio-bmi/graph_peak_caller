import logging
import numpy as np


class SparseValues:
    def __init__(self, indices, values, sanitize=False):
        self.indices = np.asanyarray(indices)
        self.values = np.asanyarray(values)
        if sanitize:
            self._sanitize()
        self.track_size = None

    def _sanitize(self):
        unique_mask = np.r_[self.indices[1:] != self.indices[:-1], True]
        self.indices = self.indices[unique_mask]
        self.values = self.values[unique_mask]

        values_mask = np.r_[True, self.values[1:] != self.values[:-1]]
        # values_mask = np.ediff1d(self.values, to_begin=1) != 0
        self.indices = self.indices[values_mask]
        self.values = self.values[values_mask]

    @classmethod
    def from_sparse_files(cls, file_base_name):
        indices = np.load(file_base_name + "_indexes.npy")
        values = np.load(file_base_name + "_values.npy")
        size = indices[-1]
        obj = cls(indices[:-1], values)
        obj.track_size = size
        return obj

    def to_sparse_files(self, file_base_name):
        np.save(file_base_name + "_indexes.npy",
                np.r_[self.indices, self.track_size])
        np.save(file_base_name + "_values.npy", self.values)
        logging.info("Wrote to %s_indexes/values.npy" % file_base_name)

    def __repr__(self):
        return "SV(%s, %s)" % (self.indices, self.values)

    def __eq__(self, other):
        if not np.all(self.indices == other.indices):
            return False

        if not np.all(self.values == other.values):
            return False

        if self.track_size != other.track_size:
            return False

        return True

    def to_bed_graph(self, filename):
        logging.warning("Not writing to %s", filename)

    def to_bed_file(self, filename):
        logging.warning("Not writing to %s", filename)

    def threshold_copy(self, cutoff):
        values = self.values >= cutoff
        new = SparseValues(self.indices, values, sanitize=True)
        new.track_size = self.track_size
        return new

    def to_dense_pileup(self, size):
        if self.values.dtype == np.bool:
            values = self.values.astype("int")
        else:
            values = self.values
        diffs = np.ediff1d(values, to_begin=values[0])
        pileup = np.zeros(size+1, dtype=values.dtype)
        indices = self.indices[:diffs.size]
        pileup[indices] = diffs
        pileup = np.cumsum(pileup[:-1])
        if self.values.dtype == np.bool:
            pileup = pileup.astype("bool")
        return pileup

    @classmethod
    def from_dense_pileup(cls, pileup):
        changes = pileup[1:] != pileup[:-1]
        indices = np.r_[0, np.flatnonzero(changes)+1]
        values = pileup[indices]
        obj = cls(indices, values)
        obj.track_size = pileup.size
        print(obj)
        return obj


class SparseDiffs:
    def __init__(self, indices, diffs, sanitize=False):
        self._indices = np.asanyarray(indices)
        self._diffs = np.asanyarray(diffs)
        if sanitize:
            self._sanitize()

    def to_dense_pileup(self, size):
        valued = self.get_sparse_values()
        return valued.to_dense_pileup(size)
        # pileup = np.zeros(size+1, dtype=self._diffs.dtype)
        # pileup[self._indices] = self._diffs
        # return np.cumsum(pileup[:-1])

    def to_sparse_files(self, file_base_name):
        np.save(file_base_name + "_indexes.npy",
                self._indices)
        np.save(file_base_name + "_values.npy", self._diffs)

    def to_bed_graph(self, filename):
        logging.warning("Not writing to %s", filename)

    def __repr__(self):
        return "SD(%s, %s)" % (self._indices, self._diffs)

    def clean(self):
        sparse_values = self.get_sparse_values()
        self._indices = sparse_values.indices
        self._diffs = np.ediff1d(sparse_values.values,
                                 to_begin=sparse_values.values[0])

    def clip_min(self, min_value):
        values = np.cumsum(self._diffs)
        np.clip(values, min_value, None, values)
        self._diffs[1:] = np.diff(values)
        self._diffs[0] = values[0]
        self._sanitize()

    def get_sparse_values(self):
        args = np.argsort(self._indices, kind="mergesort")
        diffs = self._diffs[args]
        values = np.cumsum(diffs)
        return SparseValues(self._indices[args], values, sanitize=True)


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

    @classmethod
    def from_pileup(cls, pileup, node_indexes):
        indices = np.r_[pileup.starts, pileup.ends, node_indexes]
        diffs = np.r_[np.ones(len(pileup.starts)),
                      -1*np.ones(len(pileup.ends)),
                      pileup.node_starts]
        args = np.argsort(indices)
        indices = indices[args]
        diffs = diffs[args]
        return cls(indices, diffs)

    @classmethod
    def from_dense_pileup(cls, pileup):
        diffs = np.ediff1d(pileup, to_begin=1)
        changes = np.flatnonzero(diffs)
        diffs = diffs[changes]
        diffs[0] = pileup[0]
        return cls(changes, diffs)

    def __itruediv__(self, scalar):
        self._diffs /= scalar
        return self

    def __imul__(self, scalar):
        self._diffs *= scalar
        return self

    def apply_binary_func(self, func, other, return_values=False):
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
        ret = func(values[0], values[1])[unique_mask]
        if return_values:
            return SparseValues(new_indexes, ret, sanitize=True)
        return SparseDiffs(
            new_indexes, np.ediff1d(ret, to_begin=ret[0]))


def _get_random_vec(n_points, size):
    from numpy.random import rand, randint
    return SparseDiffs(np.sort(randint(0, size, n_points)),
                       rand(n_points))

if __name__ == "__main__":
    import cProfile
    a = _get_random_vec(20000000, 500000000)
    b = _get_random_vec(20000000, 500000000)
    cProfile.run("a.maximum(b)")
