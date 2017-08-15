import sys
import pyBigWig
import numpy as np
import pybedtools
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def peaks_are_treshold_of_signal(signal_file_name, peaks_file_name):
    signal_file = pyBigWig.open(signal_file_name)
    peaks_file = pybedtools.BedTool(peaks_file_name)

    #signals = signal_file.values("chr7", 1000000, 1000100)
    #print(signals)
    #return

    for peak in peaks_file:
        chrom = peak.chrom.replace("chrom", "")
        start = peak.start
        stop = peak.stop
        #print(chrom, start, stop)

        signal_at_start = signal_file.values(chrom, start-1, start + 2)
        print(signal_at_start)


def get_long_segments(signal_file_name, peaks_file_name, min_segment_length=1e6):
    signal_file = pyBigWig.open(signal_file_name)
    peaks_file = pybedtools.BedTool(peaks_file_name)
    long_segments = []
    for peak in peaks_file:
        chrom = peak.chrom.replace("chrom", "")
        start = peak.start
        end = peak.stop
        if end - start > min_segment_length:
            signal_value = signal_file.values(chrom, start, end)[0]
            long_segments.append([end- start, signal_value])

    return long_segments




if __name__ == "__main__":
    #peaks_are_treshold_of_signal("../data/sample1_signalpvalue_1.bigwig_", "../data/sample1_peaks_1.bed")

    longest_segment("../data/sample1_signal_1.bedGraph")