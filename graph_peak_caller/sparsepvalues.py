import pickle
import numpy as np
import os
from scipy.stats import poisson
import logging
import math
from .sparsediffs import SparseValues


class PValuesFinder:
    def __init__(self, sample_pileup, control_pileup):
        self.sample = sample_pileup
        self.control = control_pileup

    def get_p_values_pileup(self):
        baseEtoTen = np.log(10)

        def clean_p_values(counts, lambdas):
            p_values = poisson.logsf(counts, lambdas)
            p_values /= -baseEtoTen
            p_values[counts == 0] = 0
            p_values[np.isinf(p_values)] = 1000
            return p_values

        p_values = self.sample.apply_binary_func(
            clean_p_values, self.control,
            return_values=True)

        return p_values


class PToQValuesMapper:

    def __init__(self, p_values, cum_counts):
        self.p_values = p_values
        self.cum_counts = cum_counts

    @classmethod
    def __read_file(cls, file_name):
        indices = np.load(file_name + "_indexes.npy")
        values = np.load(file_name + "_values.npy")
        return indices, values

    @classmethod
    def _from_subcounts(cls, p_values, counts):
        p_values = p_values.ravel()
        counts = counts.ravel()
        args = np.argsort(p_values)[::-1]
        sorted_ps = p_values[args]
        sorted_lens = counts[args]
        cum_counts = np.cumsum(sorted_lens)
        changes = np.ediff1d(sorted_ps, to_end=1) != 0
        cum_counts = cum_counts[changes]
        return cls(sorted_ps[changes], cum_counts)

    @classmethod
    def from_p_values_pileup(cls, p_values):
        logging.info("Creating mapping from p value dense pileup")
        sub_counts = cls.__get_sub_counts(p_values)
        # sub_counts = np.ediff1d(
        #     p_values.indices,
        #     to_end=p_values.track_size-p_values.indices[-1])
        return cls._from_subcounts(p_values.values, sub_counts)

    @classmethod
    def __get_sub_counts(cls, sparse_values):
        return np.ediff1d(
            sparse_values.indices,
            to_end=sparse_values.track_size-sparse_values.indices[-1])

    @classmethod
    def from_files(cls, base_file_name):
        search = base_file_name
        logging.info("Searching for files starting with %s" % search)
        files = (f for f in os.listdir()
                 if f.startswith(search) and
                 "pvalues" in f and f.endswith("_indexes.npy"))
        sub_counts = []
        p_values = []
        for filename in files:
            base_file_name = filename.replace("_indexes.npy", "")
            logging.info("Reading p values from file %s" % base_file_name)
            chr_p_values = SparseValues.from_sparse_files(base_file_name)
            sub_counts.append(cls.__get_sub_counts(chr_p_values))
            p_values.append(chr_p_values.values)
            assert sub_counts[-1].size == p_values[-1].size

        return cls._from_subcounts(
            np.concatenate(p_values),
            np.concatenate(sub_counts))

    def get_p_to_q_values(self):
        logN = np.log10(self.cum_counts[-1])
        q_values = self.p_values + np.log10(
            1+np.r_[0, self.cum_counts[:-1]])-logN
        q_values[0] = max(0, q_values[0])
        q_values = np.minimum.accumulate(q_values)
        d = dict(zip(self.p_values, q_values))
        d[0] = 0
        return d

    def to_file(self, base_name):
        with open(base_name + 'p2q.pkl', 'wb') as f:
            pickle.dump(self._p_to_q_values, f, pickle.HIGHEST_PROTOCOL)


class QValuesFinder:
    def __init__(self, p_values_pileup, p_to_q_values):
        assert isinstance(p_to_q_values, dict)
        self.p_values = p_values_pileup
        self.p_to_q_values = p_to_q_values

    def get_q_values(self):
        q_values = SparseValues(
            self.p_values.indices,
            self.get_q_array_from_p_array(self.p_values.values))
        return q_values

    def get_q_array_from_p_array(self, p_values):
        assert isinstance(p_values, np.ndarray)

        def translation(x):
            # if abs(x) < 1e-9:
            #    return 0
            # if math.isnan(x):
            #    return 0
            # x = "%.7f" % x
            # if x not in self.p_to_q_values:
            #     print(self.p_to_q_values)
            #     print(x)
            #     logging.error("P value not found in mapping dict. Could be due to rounding errors.")
            return self.p_to_q_values[x]

        trans = np.vectorize(translation, otypes=[np.float])
        new_values = trans(p_values)
        return new_values
