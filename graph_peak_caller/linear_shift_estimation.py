import numpy as np
from scipy.signal import argrelextrema
from scipy import signal
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s, %(levelname)s: %(message)s")


def moving_average(a, n=3):
    # Source: https://stackoverflow.com/a/14314054/1030104
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

class LinearShiftEstimator(object):

    def __init__(self, forward_values=None, reverse_values=None, min_pileup_value=0, window_size=200):
        self.forward_values = forward_values
        self.reverse_values = reverse_values
        self.min_pileup_value = min_pileup_value
        self.window_size = window_size

    def run(self):
        self.remove_noise()
        forward_peaks = self.find_subpeak_indices(self.forward_values)
        reverse_peaks = self.find_subpeak_indices(self.reverse_values)
        #print("Forward peaks: %s " % forward_peaks)
        #print("Reverse peaks: %s " % reverse_peaks)

        forward, reverse = self.find_subpeak_pairs(forward_peaks, reverse_peaks)

        print("N pairs: %d" % len(forward))
        #print(forward)
        #print(reverse)
        fragment_length = self.predict_fragment_length(forward, reverse)
        print("Avg fragment_length: %.2f" % fragment_length)
        return fragment_length

    def remove_noise(self):
        logging.info("Removing noise")
        self.forward_values[np.where(self.forward_values <= self.min_pileup_value)[0]] = 0
        self.reverse_values[np.where(self.reverse_values <= self.min_pileup_value)[0]] = 0

    def find_subpeak_indices(self, values):
        diffs = np.ediff1d(values, to_begin=[values[0]])
        up = list(np.where(diffs > 0)[0])
        down = list(-np.where(diffs < 0)[0])
        all_sorted = sorted(up + down, key=lambda x: abs(x))
        all_sorted = np.array(all_sorted)
        up_next_down = np.where(np.logical_and(all_sorted[:-1] > 0,all_sorted[1:] < 0))[0]
        subpeak_indices = (all_sorted[up_next_down] + np.abs(all_sorted[up_next_down+1]))/2
        return subpeak_indices

    def find_subpeak_pairs(self, forward_peaks, reverse_peaks):
        forward_peaks = list(forward_peaks)
        reverse_peaks = list(-reverse_peaks)

        all_sorted = sorted(forward_peaks + reverse_peaks, key=lambda x: abs(x))
        all_sorted = np.array(all_sorted)
        forward_and_next_reverse = np.where(np.logical_and(all_sorted[:-1] > 0, all_sorted[1:] < 0))[0]

        forward_indices = all_sorted[forward_and_next_reverse]
        reverse_indices = -all_sorted[forward_and_next_reverse+1]

        return forward_indices, reverse_indices

    def predict_fragment_length(self, forward_indices, reverse_indices):
        sizes = reverse_indices - forward_indices
        sizes = sizes[np.where(sizes > 20)]
        sizes = sizes[np.where(sizes < 200)]
        #for size in sizes:
        #    print(size)

        return np.median(sizes)


if __name__ == "__main__":

    def benchmark():
        forward_values = np.zeros(200000000)
        reverse_values = np.zeros(200000000)

        random_starts = np.random.randint(0, 20000000, size=100000, dtype=np.int64)
        for start in random_starts:
            forward_values[start:start + 50] += 1

        random_starts = np.random.randint(0, 20000000, size=100000, dtype=np.int64)
        for start in random_starts:
            reverse_values[start:start + 50] += 1

        estimator = LinearShiftEstimator(forward_values, reverse_values)
        estimator.run()

    def test_real_data():
        forward = np.zeros(300000000)
        reverse = np.zeros(300000000)
        f = open("../tests/ctcf_reads_chr1.bed")
        i = 0
        for line in f:
            if i % 100000 == 0:
                logging.info("Read %d" % i)

            i += 1

            l = line.split()
            start = int(l[1])
            end = int(l[2])
            strand = l[5]

            if strand == "+":
                forward[start:end] += 1
            else:
                reverse[start:end] += 1


        estimator = LinearShiftEstimator(forward, reverse, min_pileup_value=5)
        estimator.run()

    test_real_data()