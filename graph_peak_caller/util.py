import sys
import pyBigWig
import numpy as np
import pybedtools
from pybedtools import BedTool
from offsetbasedgraph import IntervalCollection
from offsetbasedgraph.graphtraverser import GraphTraverserUsingSequence
import offsetbasedgraph as obg
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

def bed_intervals_to_graph(obg_graph, linear_path_interval, bed_file_name, graph_start_offset):
    peaks = BedTool(bed_file_name)
    intervals_on_graph = []
    for peak in peaks:
        start = peak.start - graph_start_offset
        end = peak.end - graph_start_offset
        intervals_on_graph.append(linear_path_interval.get_subinterval(start, end))

    return intervals_on_graph

def fasta_sequence_to_linear_path_through_graph(linear_sequence_fasta_file, sequence_retriever, ob_graph, start_node):
    search_sequence = open(linear_sequence_fasta_file).read()
    print("Length of search sequence: %d" % len(search_sequence))
    traverser = GraphTraverserUsingSequence(ob_graph, search_sequence, sequence_retriever)
    traverser.search_from_node(start_node)
    linear_path_interval = traverser.get_interval_found()
    #IntervalCollection([linear_path_interval]).to_file("linear_path", text_file=True)
    return linear_path_interval


def get_linear_paths_in_graph(ob_graph, vg_graph, write_to_file_name=None):
    intervals = []
    for path in vg_graph.paths:
        #print(path.__dict__.keys())
        print(path.name)
        obg_interval = path.to_obg(ob_graph=ob_graph)
        print(obg_interval.length())
        obg_interval.name = path.name
        intervals.append(obg_interval)

    if write_to_file_name is not None:
        collection = obg.IntervalCollection(intervals)
        collection.to_file(write_to_file_name, text_file=True)

    return intervals


def sparse_maximum(indices1, values1, indices2, values2, genome_size):

    if indices1[0] != 0:
        np.insert(indices1, 0, 0)
        np.insert(values1, 0, 0)

    if indices2[0] != 0:
        np.insert(indices2, 0, 0)
        np.insert(values2, 0, 0)

    indices1 = np.insert(indices1, len(indices1), genome_size)
    indices2 = np.insert(indices2, len(indices2), genome_size)


    a = indices1[:-1]*2
    b = indices2[:-1]*2+1
    all_idxs = np.concatenate([a, b])
    all_idxs.sort()
    vi_list = [values1, values2]
    values_list = []
    for i, vals in enumerate(vi_list):
        idxs = np.nonzero((all_idxs % 2) == i)[0]
        all_values = vals
        #print("Idx: %s" % idxs)
        #print("Values %i: %s" % (i, all_values))
        value_diffs = np.diff(all_values)
        #print("Diffs: %s" % value_diffs)
        values = np.zeros(all_idxs.shape)
        values[idxs[1:]] = value_diffs
        values[idxs[0]] = all_values[0]
        #print("  values after diffs: %s" % values)
        values_list.append(values.cumsum())
    #print("Values list:")
    #print(values_list)
    values = np.maximum(values_list[0], values_list[1])
    idxs = all_idxs // 2
    empty_ends = np.nonzero(np.diff(idxs) == 0)[0]
    max_values = np.maximum(values[empty_ends], values[empty_ends+1])
    values[empty_ends+1] = max_values
    values[empty_ends] = max_values

    indices = idxs
    values = values

    indices, values = sanitize_indices_and_values(indices, values)
    return indices, values

def sanitize_indices_and_values(indices, values):

    new_indices = []
    new_values = []

    prev = None
    for i, index in enumerate(indices):
        if values[i] != prev:
            new_indices.append(index)
            new_values.append(values[i])
        prev = values[i]

    return new_indices, new_values

if __name__ == "__main__":
    values = get_average_signal_values_within_peaks("../data/sample1_signal_1.bigwig", "../data/sample1_peaks_1.bed")
    values = values[np.logical_not(np.isnan(values))]
    values = values[values < 8]
    print("Mean", np.mean(values))
    print("Std:", np.std(values))
    plt.hist(values, 100)
    plt.show()
    #longest_segment("../data/sample1_signal_1.bedGraph")
