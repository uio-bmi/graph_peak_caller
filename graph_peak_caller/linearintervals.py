import numpy as np


class LinearIntervalCollection(object):

    def __init__(self, starts, ends):
        self.starts = starts
        self.ends = ends

    def extend(self, extension_size):
        extended_starts = (self.starts + self.ends)/2 - extension_size
        extended_ends = (self.starts + self.ends)/2 + extension_size
        return self.__class__(extended_starts, extended_ends)

    def n_basepairs_covered(self):
        return np.sum(self.ends - self.starts)