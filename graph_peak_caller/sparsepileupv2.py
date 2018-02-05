import logging
from itertools import chain
import numpy as np
from scipy.stats import poisson
from collections import defaultdict
from .pileup import Pileup
from .pileupcleaner2 import PeaksCleaner, HolesCleaner
from .subgraphcollection import SubgraphCollection
from .eventsorter import DiscreteEventSorter
from offsetbasedgraph import Interval, IntervalCollection
import pickle
from .sparsepileup import SparseAreasDict, starts_and_ends_to_sparse_pileup, intervals_to_start_and_ends
from memory_profiler import profile

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

