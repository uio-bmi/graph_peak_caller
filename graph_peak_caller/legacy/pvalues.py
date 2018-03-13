import pickle
import numpy as np
import os
from collections import Counter
from .densepileup import DensePileup
from scipy.stats import poisson
import logging
import math


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

        # Set p-value to 0 where sample is 0
        zero_indices = np.where(self.sample.data._values == 0)[0]
        p_values_array[zero_indices] = 0

        p_values_pileup = DensePileup(self.sample.graph)
        p_values_pileup.set_new_values(p_values_array)
        p_values_pileup.data._touched_nodes = self.sample.data._touched_nodes
        return p_values_pileup


class PToQValuesMapper:

    def __init__(self, counter):
        self._counter = counter

    @classmethod
    def __read_file(cls, file_name):
        indices = np.load(file_name + "_indexes.npy")
        values = np.load(file_name + "_values.npy")
        return indices, values

    @classmethod
    def from_p_values_dense_pileup(cls, p_values):
        logging.info("Creating mapping from p value dense pileup")
        values = p_values.data._values
        sorted_p_values = sorted(values, reverse=True)
        unique, counts = np.unique(sorted_p_values, return_counts=True)
        counter_dict = {unique[i]: counts[i]
                        for i in range(0, len(counts))}
        return cls(counter_dict)

    @classmethod
    def from_files(cls, base_file_name):
        search = base_file_name
        logging.info("Searching for files starting with %s" % search)
        files = (f for f in os.listdir()
                 if f.startswith(search) and "pvalues" in f and f.endswith("_indexes.npy"))
        counter = Counter()
        for filename in files:
            base_file_name = filename.replace("_indexes.npy", "")  # Get file name base
            logging.info("Reading p values from file %s" % base_file_name)
            indices, values = cls.__read_file(base_file_name)
            counts = np.diff(indices)
            for key, val in zip(values, counts):
                counter[key] += val
            #counter.update(dict(zip(values, counts)))
        return cls(counter)

    def get_p_to_q_values(self):
        logging.info("Mapping p to q values")
        p_to_q_values = {}
        rank = 1
        logN = np.log10(sum(self._counter.values()))
        pre_q = None
        i = 0
        for p_value in reversed(sorted(self._counter.keys())):
            if i % 5000 == 0:
                logging.info("Done mapping p value %d/%d" % (i, len(self._counter.keys())))
            i += 1

            value_count = self._counter[p_value]
            q_value = p_value + (np.log10(rank) - logN)
            if rank == 1:
                q_value = max(0.0, q_value)
            else:
                q_value = max(0.0, min(pre_q, q_value))
            p_to_q_values["%.7f" % p_value] = q_value
            pre_q = q_value
            rank += value_count

        self.p_to_q_values = p_to_q_values
        logging.info("Done creating mapping")
        return p_to_q_values

    def to_file(self, base_name):
        with open(base_name + 'p2q.pkl', 'wb') as f:
            pickle.dump(self._p_to_q_values, f, pickle.HIGHEST_PROTOCOL)


class QValuesFinder:
    def __init__(self, p_values_pileup, p_to_q_values):
        assert isinstance(p_to_q_values, dict)
        self.p_values = p_values_pileup
        self.p_to_q_values = p_to_q_values

    def get_q_values(self):
        q_values_pileup = DensePileup(self.p_values.graph)
        new_values = self.get_q_array_from_p_array(
                        self.p_values.data._values)
        q_values_pileup.set_new_values(new_values)
        q_values_pileup.data._touched_nodes = self.p_values.data._touched_nodes
        return q_values_pileup

    def get_q_array_from_p_array(self, p_values):
        assert isinstance(p_values, np.ndarray)

        def translation(x):
            if abs(x) < 1e-9:
                return 0
            if math.isnan(x):
                return 0
            x = "%.7f" % x
            if x not in self.p_to_q_values:
                print(self.p_to_q_values)
                logging.error("P value not found in mapping dict. Could be due to rounding errors.")
            return self.p_to_q_values[x]

        trans = np.vectorize(translation, otypes=[np.float])
        new_values = trans(p_values)
        return new_values
