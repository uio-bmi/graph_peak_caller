import numpy as np


class LinearIntervalCollection(object):

    def __init__(self, starts, ends):
        self.starts = np.asanyarray(starts)
        self.ends = np.asanyarray(ends)
        self.n_intervals = self.starts.size

    def __str__(self):
        return str(self.starts) + "\n" + str(self.ends)

    def extend_mid(self, extension_size):
        extended_starts = (self.starts + self.ends)/2 - extension_size
        extended_ends = (self.starts + self.ends)/2 + extension_size
        return self.__class__(extended_starts, extended_ends)

    def extend(self, extension_size):
        extended_starts = self.starts - extension_size
        extended_ends = self.starts + extension_size
        return self.__class__(extended_starts, extended_ends)

    def extend_np(self, extension_size):
        return np.add.outer(np.array([-extension_size, extension_size]),
                            self.starts)

    def n_basepairs_covered(self):
        return np.sum(np.abs(self.ends - self.starts))
