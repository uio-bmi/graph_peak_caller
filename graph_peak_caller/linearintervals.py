import numpy as np


class LinearIntervalCollection(object):

    def __init__(self, starts, ends):
        self.starts = np.array(starts)
        self.ends = np.array(ends)
        self.n_intervals = self.starts.size

    def __str__(self):
        return str(self.starts) + "\n" + str(self.ends)

    def extend(self, extension_size):
        extended_starts = (self.starts + self.ends)/2 - extension_size
        extended_ends = (self.starts + self.ends)/2 + extension_size
        return self.__class__(extended_starts, extended_ends)

    def n_basepairs_covered(self):

        lengths = self.ends - self.starts
        for l in lengths:
            print(l)

        return np.sum(np.abs(self.ends - self.starts))
