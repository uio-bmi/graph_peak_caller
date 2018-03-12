import numpy as np


class LinearIntervalCollection(object):

    def __init__(self, starts, ends):
        self.starts = np.asanyarray(starts)
        self.ends = np.asanyarray(ends)
        self.n_intervals = self.starts.size

    def __eq__(self, other):
        if not np.allclose(self.starts, other.starts):
            return False
        print(self.ends)
        print(other.ends)
        return np.allclose(self.ends, other.ends)

    def __repr__(self):
        s = " ".join(str(i) for i in self.starts)
        e = " ".join(str(i) for i in self.ends)
        return "(%s:%s)" % (s, e)

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
