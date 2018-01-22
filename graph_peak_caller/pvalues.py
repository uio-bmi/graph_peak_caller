import numpy as np
from scipy.stats import poisson
import logging


class PValuesFinder:
    def __init__(self, sample_pileup, control_pileup):
        self.sample = sample_pileup
        self.control = control_pileup

    def get_p_values_array(self):
        p_values = poisson.logsf(
                                self.sample.data._values,
                                self.control.data._values)
        baseEtoTen = np.log(10)
        self.p_values_array -p_values / baseEtoTen

    def to_sparse_file(self, truncate_below=0.05, file_base_name="p_values"):
        self.p_values_array[np.where(self.p_values_array < truncate_below)] = 0
        indexes = np.where(np.ediff1d(self.p_values_array, to_begin=[0]) != 0)
        values = self.p_values_array[indexes]

        np.savetxt(file_base_name + "_indexes.npy", indexes)
        np.savetxt(file_base_name + "_values.npy", values)

        logging.info("Saved p values indexes/values to files")


class PToQValuesMapper:

    @classmethod
    def from_files(cls, base_file_name):
        pass

    def to_file(self):
        pass


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
