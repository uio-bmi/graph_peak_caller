import sys
import pyBigWig
import numpy as np
import pybedtools

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def get_average_signal_values_within_peaks(signal_file_name, peaks_bed_file_name):
    signal_file = pyBigWig.open(signal_file_name)
    peaks_file = pybedtools.BedTool(peaks_bed_file_name)
    long_segments = []
    max_segments = 200000
    i = 0
    for peak in peaks_file:
        chrom = peak.chrom.replace("chrom", "")
        start = peak.start
        end = peak.stop
        signal_values = signal_file.values(chrom, start, end)
        #signal_value = np.log(np.mean(signal_values))
        signal_value = np.mean(signal_values)
        long_segments.append(signal_value)
        #print("Found %d-%d with value %.5f" % (start, end, signal_value))
        i += 1
        if i > max_segments:
            break

    return np.array(long_segments)


def longest_segment(file_name):
    longest = 0
    for line in open(file_name):
        l = line.split()
        start = int(l[1])
        end = int(l[2])
        #print(end - start)
        if end - start > longest:
            print("Found longer %d at %d-%d with value %.6f" % (end-start, start, end, float(l[3])))
            longest = end - start

    print(longest)


if __name__ == "__main__":
    values = get_average_signal_values_within_peaks("../data/sample1_signal_1.bigwig", "../data/sample1_peaks_1.bed")
    values = values[np.logical_not(np.isnan(values))]
    values = values[values < 8]
    print("Mean", np.mean(values))
    print("Std:", np.std(values))
    plt.hist(values, 100)
    plt.show()
    #longest_segment("../data/sample1_signal_1.bedGraph")