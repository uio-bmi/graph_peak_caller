import sys
import pyBigWig
import numpy as np
import matplotlib.pyplot as plt

CHUNK_SIZE = 100000
SUBSAMPLING = 10


class BigWigHistogram(object):

    def __init__(self, file_name, nbins=1000, transform=np.log):
        self.transform = transform
        self.file_name = file_name
        self.file_object = pyBigWig.open(file_name)
        self.chromosomes = self.file_object.chroms()
        self.histogram = np.zeros(nbins, dtype="int64")
        self.nbins = nbins
        header = self.file_object.header()
        print(header)
        self.max_signal = self.transform(header["maxVal"])
        self.min_signal = self.transform(header["minVal"])
        self.range = (self.min_signal, (self.min_signal+self.max_signal)/2)
        self.bin_size = (self.max_signal-self.min_signal)/self.nbins
        print(self.range)

    def create_histogram(self):
        for chromosome in self.chromosomes:
            self._get_chromosome_histogram(chromosome)

    def _update_histogram(self, chunk):
        t_chunk = self.transform(chunk)
        local_histogram = np.histogram(
            t_chunk, self.nbins, self.range)[0]
        self.histogram += local_histogram

    def _get_chromosome_histogram(self, chromosome):
        length = self.chromosomes[chromosome]
        n_chunks = length//CHUNK_SIZE+1
        for i in range(n_chunks):
            if i % SUBSAMPLING:
                continue
            start = i*CHUNK_SIZE
            end = min(length, start+CHUNK_SIZE)
            print("\t Reading chunk %s of %s (%s:%s)" % (i, n_chunks, start, end))
            chunk = np.array(self.file_object.values(
                chromosome, start, end))
            if any(np.isnan(chunk)):
                print(sum(np.isnan(chunk)))
                print (chromosome, end, length)
                continue
            self._update_histogram(chunk)

    def save(self):
        new_file_name = self.file_name.replace("bigwig", "npy").replace("bigWig", "npy")
        assert new_file_name != self.file_name
        np.save(new_file_name, self.histogram)

    def plot(self):
        plt.plot(self.histogram)
        plt.show()


class IntervalHistogram(BigWigHistogram):
    def _get_chromosome_histogram(self, chromosome):
        print("#", chromosome)
        length = self.chromosomes[chromosome]
        intervals = self.file_object.intervals(chromosome, 0, length)
        for interval in intervals:
            interval_length = interval[1]-interval[0]
            value = self.transform(interval[2])
            bin_number = int((value-self.range[0])/self.bin_size)
            if bin_number > self.nbins-1:
                print(bin_number, value)
                bin_number = self.nbins-1
            self.histogram[bin_number] += interval_length

    def plot_all_lengths(self):
        for chromosome in self.chromosomes:
            if "_" in chromosome:
                continue
            self._plot_lengths(chromosome)

    def _plot_lengths(self, chromosome):
        print("#", chromosome)
        length = self.chromosomes[chromosome]
        intervals = self.file_object.intervals(chromosome, 0, length)
        intervals = [i for i in intervals if i[2] != 0]
        lengths = np.zeros(len(intervals))
        values = np.zeros(len(intervals))
        for i, interval in enumerate(intervals):
            lengths[i] = log(interval[1]-interval[0])
            values[i] = interval[2]  # /0.4621149897575
            if lengths[i] > 100:
                print(values[i], lengths[i])
        plt.scatter(values, lengths)
        plt.show()


class BedGraphHistogram():
    def __init__(self, file_name):
        self.file_object = open(file_name)

    def _parse_line(self, line):
        parts = line.split()
        length = int(parts[2]) - int(parts[1])
        value = float(parts[3])
        return (value, length)

    def plot_all_lengths(self):
        pairs = [self._parse_line(line) for line in
                 self.file_object.readlines() if not line.startswith("#")]
        values = np.array([pair[0] for pair in pairs])
        lengths = np.array([pair[1] for pair in pairs])
        plt.scatter(values, lengths)
        plt.show()

if __name__ == "__main__":
    np.seterr("raise")
    file_name = sys.argv[1]
    if file_name.lower().endswith("bigwig"):
        histogram = IntervalHistogram(file_name,
                                      transform=lambda x: np.log(x+1))
    elif file_name.lower.endswith("bdg"):
        histogram = BedGraphHistogram(file_name)
    histogram.plot_all_lengths()
