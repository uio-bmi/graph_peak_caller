import numpy as np


class SimpleValuedIndexes():
    def __init__(self, indexes, values):
        self.indexes = indexes
        self.values = values

    def sum(self):
        lengths = np.diff(self.indexes)
        return np.sum(lengths*self.values)


class RpScore:
    def __init__(self, max_score, sum_score):
        self.sum_score = sum_score
        self.max_score = max_score

    def __getitem__(self, item):
        if item == 0:
            return self.max_score
        elif item == 1:
            return self.sum_score
        else:
            raise NotImplementedError()

    def sum(self):
        return self.sum_score

    def max_value(self):
        return self.max_score

    @classmethod
    def from_valued_indexes(cls, vi):
        return cls(np.max(vi.all_values()), vi.sum())

