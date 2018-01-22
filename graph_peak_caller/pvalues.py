import pickle
import numpy as np
import os
from collections import Counter
from .densepileup import DensePileup
from scipy.stats import poisson
import logging


class PValuesFinder:
    def __init__(self, sample_pileup, control_pileup):
        self.sample = sample_pileup
        self.control = control_pileup

    def get_p_values_pileup(self):
        p_values = poisson.logsf(
                                self.sample.data._values,
                                self.control.data._values)
        baseEtoTen = np.log(10)
        p_values_array = -p_values / baseEtoTen
        p_values_pileup = DensePileup(self.sample.graph)
        p_values_pileup.set_new_values(p_values_pileup)


class PToQValuesMapper:

    def __init__(self, counter, base_name):
        self._counter = counter
        self._base_name = base_name

    @classmethod
    def __read_file(cls, file_name):
        indices = np.load(file_name + "_indices.npy")
        values = np.load(file_name + "_values.npy")
        return indices, values

    @classmethod
    def from_files(cls, base_file_name):
        files = (f for f in os.listdir()
                 if f.startswith(base_file_name + "_chr"))
        counter = Counter()
        for filename in files:
            indices, values = cls.__read_file(filename)
            counts = np.diff(indices)
            counter.update(dict(zip(values, counts)))
        return cls(counter, base_file_name)

    def get_p_to_q_values(self):
        p_to_q_values = {}
        rank = 1
        logN = np.log10(sum(self._counter.values()))
        pre_q = None
        for p_value in reversed(sorted(self._counter.keys())):
            value_count = self.counter[p_value]
            q_value = p_value + (np.log10(rank) - logN)
            if rank == 1:
                q_value = max(0.0, q_value)
            else:
                q_value = max(0.0, min(pre_q, q_value))

            p_to_q_values[p_value] = q_value
            pre_q = q_value
            rank += value_count

        self.p_to_q_values = p_to_q_values

    def to_file(self):
        with open(self.base_name + 'p2q.pkl', 'wb') as f:
            pickle.dump(self._counter, f, pickle.HIGHEST_PROTOCOL)


class QValuesFinder:
    def __init__(self, p_values_pileup, p_to_q_values):
        assert isinstance(p_to_q_values, dict)
        self.p_values = p_values_pileup
        self.p_to_q_values = p_to_q_values

    def get_q_values(self):
        new_values = self.get_q_array_from_p_array(
                        self.p_values.data._values)
        self.p_values.data._values = new_values
        return self.p_values

    def get_q_array_from_p_array(self, p_values):
        assert isinstance(p_values, np.ndarray)

        def translation(x):
            if x == 0:
                return 0
            return self.p_to_q_values[x]

        trans = np.vectorize(translation, otypes=[np.float])
        new_values = trans(self.p_values)
        return new_values
