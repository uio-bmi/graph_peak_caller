from collections import defaultdict
import matplotlib.pyplot as plt
import pyvg
from pyvg.sequences import SequenceRetriever
import offsetbasedgraph as obg
from offsetbasedgraph.graphtraverser import GraphTraverserUsingSequence
from graph_peak_caller.peakcollection import PeakCollection
from offsetbasedgraph import IntervalCollection, DirectedInterval
from graph_peak_caller.subgraphcollection import SubgraphCollection
from graph_peak_caller.util import LinearRegion, get_linear_paths_in_graph
from graph_peak_caller.sparsepileup import SparsePileup, ValuedIndexes
from peakscomparer import PeaksComparer
from compareMacs import ValuedInterval
import numpy as np

MHC_REGION = LinearRegion("chr6", 28510119, 33480577)


def create_linear_path(ob_graph, vg_graph):
    linear_paths = get_linear_paths_in_graph(ob_graph, vg_graph, "linear_maps")
    ref_path = linear_paths["ref"].to_indexed_interval()
    return ref_path


def create_linear_peaks_from_bed(linear_sequence_fasta_file, peaks_bed_file,
                                 obg_graph_file_name, vg_graph_file_name,
                                 start_node,
                                 region):

    ob_graph = obg.GraphWithReversals.from_file(obg_graph_file_name)
    search_sequence = open(linear_sequence_fasta_file).read()
    sequence_retriever = SequenceRetriever.from_vg_graph(vg_graph_file_name)
    traverser = GraphTraverserUsingSequence(
        ob_graph, search_sequence, sequence_retriever)
    traverser.search_from_node(start_node)
    linear_path_interval = traverser.get_interval_found()
    IntervalCollection([linear_path_interval]).to_file("linear_path.intervalcollection",
                                                       text_file=True)
    print("Length")
    print(linear_path_interval.length())
    print(linear_path_interval.region_paths[0])
    print(linear_path_interval.start_position)
    print(linear_path_interval.end_position)

    linear_peaks = PeakCollection.create_from_linear_intervals_in_bed_file(
                        obg_graph_file_name,
                        linear_path_interval,
                        peaks_bed_file,
                        region.start,
                        region.end)

    linear_peaks.to_file("linear_peaks.intervalcollection", text_file=True)


class SubgraphAnalyser(object):

    def __init__(self, graph, subgraphs_file_name):
        self.graph = graph
        self.subgraphs = SubgraphCollection.from_pickle(
            subgraphs_file_name, graph=graph)

    def print_sizes(self):
        sizes = []
        for subgraph in self.subgraphs:
            print(subgraph.n_basepairs())
            sizes.append(subgraph.n_basepairs())

        plt.hist(sizes, bins=100)
        plt.show()


class SubgraphComparer():
    # Compare subgraph vs linear peaks
    def __init__(self, graph, subgraphs_file_name, peaks_file_name):
        self.graph = graph
        IntervalCollection.interval_class = DirectedInterval
        self.subgraphs = SubgraphCollection.from_pickle(
            subgraphs_file_name, graph=graph)
        self.peaks = PeakCollection.create_list_from_file(
            peaks_file_name, graph=graph)

    def check_peaks_in_subgraphs(self):
        n_in_subgraphs = 0
        for peak in self.peaks:
            print(peak)
            if self.subgraphs.contains_interval(peak):
                n_in_subgraphs += 1

        print("%d peaks are in subgraphs" % n_in_subgraphs)


class AlignmentsAnalyser(object):
    def __init__(self, vg_graph, vg_gam_file_name, ob_graph,
                 linear_path_interval_file_name):
        self.graph = ob_graph
        self.vg_graph = vg_graph
        print("Reading reads")
        self.reads = vg_gam_file_to_interval_list(
            vg_graph, vg_gam_file_name, ob_graph, max_intervals=10000)
        print("Number of reads: %d" % len(self.reads))

        self.linear_path = IntervalCollection.create_list_from_file(
            linear_path_interval_file_name, self.graph).intervals[0]

    def count_alignments_on_linear_path(self):
        n_on_path = 0
        n_intersects = 0
        n_not_on_path = 0
        i = 0
        for read in self.reads:
            if i % 100 == 0:
                print("#%d" % i)
            i += 1
            if self.linear_path.contains_in_order_any_direction(read):
                n_on_path += 1
            elif self.linear_path.intersects(read) or self.linear_path.intersects(read.get_reverse()):
                n_intersects += 1
            else:
                n_not_on_path += 1

        print("N on path: %d" % n_on_path)
        print("N intersects: %d" % n_intersects)
        print("N not on path: %d" % n_not_on_path)

    def count_alignments_on_linear_paths(self):
        paths = get_linear_paths_in_graph(self.graph, self.vg_graph)
        hits = defaultdict(int)

        i = 0
        for read in self.reads:

            for path in paths:
                if path.contains(read):
                    hits[path.name] += 1

            if i % 100 == 0:
                print("#%d" % i)
            i += 1

        print("N reads: %d" % i)
        for name, hit in hits.items():
            print("%s: %d " % (name, hit))


def macs_pileup_to_graph_pileup(ob_graph, linear_path, pileup_file, region):

    graph_intervals = []
    values = []
    graph_pileup = SparsePileup(ob_graph)
    for line in open(pileup_file):
        if i  % 100 == 0:
            print("Pileup line %d" % i)
        l = line.split()
        chrom = l[0]
        if chrom != region.chromosome:
            continue

        start = int(l[1])
        end = int(l[2])

        if start < region.start or end > region.end:
            continue

        value = float(l[3])

        graph_start = start - region.start
        graph_end = end - region.start

        graph_interval = linear_path.get_subinterval(graph_start, graph_end)
        graph_intervals.append(graph_interval)
        values.append(value)

    graph_pileup.set_sorted_interval_values(graph_interval, values)

    return graph_pileup
    #return SparsePileup.from_intervals(ob_graph, graph_intervals)

def check_macs_pileup_values_for_graph_peaks(graph_peaks_file_name,
                                             macs_pileup_file_Name):
    ob_graph = obg.GraphWithReversals.from_file("obgraph")
    vg_graph = pyvg.vg.Graph.create_from_file("haplo1kg50-mhc.json")
    path = create_linear_path(ob_graph, vg_graph)


def find_missing_graph_peaks():
    ob_graph = obg.GraphWithReversals.from_file("obgraph")
    vg_graph = pyvg.vg.Graph.create_from_file("haplo1kg50-mhc.json")
    path = create_linear_path(ob_graph, vg_graph)
    comparer = PeaksComparer.create_from_graph_peaks_and_linear_peaks(
        "ctcf05_peaks.narrowPeak",
        "laststepmax_paths.intervalcollection",
        ob_graph,
        path,
        MHC_REGION)
    comparer.check_similarity()


def analyse_without_control():
    ob_graph = obg.GraphWithReversals.from_file("haplo1kg50-mhc.obg")
    vg_graph = pyvg.vg.Graph.create_from_file("haplo1kg50-mhc.json")
    path = create_linear_path(ob_graph, vg_graph)
    comparer = PeaksComparer.create_from_graph_peaks_and_linear_peaks(
        "macs_without_control_peaks.narrowPeak",
        "ctcf_q50_without_control_max_paths.intervalcollection",
        ob_graph,
        path,
        MHC_REGION)
    comparer.check_similarity()


def analyse_pileups():
    ob_graph = obg.GraphWithReversals.from_file("graph.obg")
    vg_graph = pyvg.vg.Graph.create_from_file("haplo1kg50-mhc.json")
    path = create_linear_path(ob_graph, vg_graph)

    macs_pileup = macs_pileup_to_graph_pileup(ob_graph, path, "macs_without_control_treat_pileup.bdg", MHC_REGION)
    macs_pileup.to_bed_graph("macs_sample_on_graph.bdg")
    return

    #interval = obg.Interval.from_file_line('{"region_paths": [461954, 461955, 461957, 461958, 461960, 461962, 461963, 461965, 461966, 461968, 461969, 461971, 461972, 461974, 461975, 461977], "end": 64, "start": 35, "average_q_value": 86.40319972761958, "direction": 1}')
    interval = obg.Interval.from_file_line('{"region_paths": [63746, 63747, 63749, 63750, 63752, 63753, 63755, 63756, 63758, 63759, 63761, 63763, 63764, 63765, 63767, 63768, 63770, 63772, 63773, 63775, 63776, 63778, 63779, 63781, 63783, 63784, 63785, 63787, 63788, 63790, 63792, 63794, 500765, 63796, 63798, 63800, 63801, 63802, 63804, 63806, 500768, 63808, 63810, 63811, 63813, 63814, 63816, 63818, 63820, 500772, 63822, 63824], "average_q_value": 66.00805308978333, "end": 2, "start": 82, "direction": 1}')

    linear_start = path.get_offset_at_position(interval.start_position) + MHC_REGION.start
    linear_end = path.get_offset_at_position(interval.end_position) + MHC_REGION.start

    print(linear_start)
    print(linear_end)
    return

    sample = SparsePileup.from_bed_graph(ob_graph, "ctcf_q50_without_control_sample_track.bdg")
    sample2 = SparsePileup.from_bed_graph(ob_graph, "ctcf_q50_with_control_sample_track.bdg")
    control1 = SparsePileup.from_bed_graph(ob_graph, "ctcf_q50_without_control_scaled_control.bdg")
    control2 = SparsePileup.from_bed_graph(ob_graph, "ctcf_q50_with_control_scaled_control.bdg")
    print(" == Without control == ")
    for rp in interval.region_paths:
        print("Sample without:   %d: %s" % (rp, sample.data[rp]))
        print("Sample with       %d: %s" % (rp, sample2.data[rp]))
        print("Control without   %d: %s" % (rp, control1.data[rp]))
        print("Control with      %d: %s" % (rp, control2.data[rp]))
        assert sample.data[rp] == sample2.data[rp], "\n%s != \n%s" % (sample.data[rp], sample2.data[rp])

#comparer.compare_q_values_for_similar_peaks()


#comparer.plot_peak_lengths()
#comparer.check_similarity()
#comparer.check_overlap_with_linear_path()
#peaks = comparer.get_peaks_not_on_linear_path()
#comparer.peaks_to_fasta(peaks)




#comparer = SubgraphComparer(ob_graph, "real_data_peaks_as_subgraphs.pickle", "linear_peaks.intervalcollection")
#comparer.check_peaks_in_subgraphs()

#analyser = SubgraphAnalyser(ob_graph, "real_data_peaks_as_subgraphs.pickle")
#analyser.print_sizes()

#peaks = comparer.graph_peaks_on_main_path_not_in_linear()
#comparer.peaks_to_fasta(peaks, "alone_linear.peaks")


#analyser = AlignmentsAnalyser(vg_graph, "ENCFF001HNI_filtered_q60.gam", ob_graph, "linear_path.intervalcollection")  # sample reads
#analyser = AlignmentsAnalyser(vg_graph, "ENCFF001HNS_filtered_q60.gam", ob_graph, "linear_path.intervalcollection")  # Control reads
#analyser.count_alignments_on_linear_paths()
#analyser.count_alignments_on_linear_path()

if __name__ == "__main__":
    analyse_pileups()
    #find_missing_graph_peaks()
    exit()
    ob_graph = obg.GraphWithReversals.from_file("graph.obg")
    vg_graph = pyvg.vg.Graph.create_from_file("haplo1kg50-mhc.json")
    get_linear_paths_in_graph(ob_graph, vg_graph, "linear_maps")
