import numpy as np
from scipy.stats import poisson


class PValuesFinder():
    def __init__(self, sample_pileup, control_pileup):
        self.sample = sample_pileup
        self.control = control_pileup

    def get_p_values(self):
        self.p_values = poisson.logsf(self.sample.data._values,
                                self.control.data._values)
        baseEtoTen = np.log(10)
        self.p_values = -self.p_values / baseEtoTen


class PToQValuesMapper:

    @classmethod
    def from_files(cls, base_file_name):
        pass

    def to_file(self):
        pass



def get_q_from_p_values(self, p_values, p_to_q_values):
    assert isinstance(p_values, np.ndarray)
    assert isinstance(p_to_q_values, dict)

    def translation(x):
        if x == 0:
            return 0
        return self.p_to_q_values[x]

    trans = np.vectorize(translation, otypes=[np.float])
    new_values = trans(self.p_values)  # np.apply_along_axis(translation, 0, self.p_values)

    return new_values