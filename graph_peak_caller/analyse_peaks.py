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
from graph_peak_caller.sparsepileup import SparsePileup
import numpy as np

MHC_REGION = LinearRegion("chr6", 28510119, 33480577)


def create_linear_path(ob_graph, vg_graph):
    assert ob_graph is not None
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
        #"ctcf_q50_without_control_max_paths.intervalcollection",
        "ctcf_r1_max_paths.intervalcollection",
        ob_graph,
        path,
        MHC_REGION)
    comparer.check_similarity(analyse_first_n_peaks=50)


def get_mappings_overlapping_with_interval(mappings, intervals):
    i = 0
    n = defaultdict(int)
    for i, mapping in enumerate(mappings):

        for j, interval in enumerate(intervals):
            if interval.contains_position(mapping.start_position):
                n[j] += 1

    return n

def analyse_pileups_on_peaks(ob_graph, pileups_file_names, peak_intervals_file_name):
    print("Analysing peaks")
    pileups = {name: SparsePileup.from_bed_graph(ob_graph, pileup) for name, pileup in pileups_file_names.items()}
    peaks = IntervalCollection.from_file(peak_intervals_file_name, text_file=True)

    for peak in peaks:
        print()
        print("Peak %s" % peak)
        rp = peak.region_paths[0]
        for name, pileup in pileups.items():
            pileup_sum = sum(pileup.data[rp].sum() for rp in peak.region_paths)
            print("Pileup %s: %d" % (name, pileup_sum))


def vg_alignments_to_linear():
    ob_graph = obg.GraphWithReversals.from_file("haplo1kg50-mhc.obg")
    vg_graph = pyvg.vg.Graph.create_from_file("haplo1kg50-mhc.json")
    path = create_linear_path(ob_graph, vg_graph)
    analyser = AlignmentsAnalyser(vg_graph, "ENCFF001HNI_haplo1kg50-mhc_filtered_q50.gam", ob_graph, path)  # sample reads
    #linear = analyser.to_linear_alignments()
    #collection = IntervalCollection(linear)
    #collection.to_file("graph_reads_on_linear2.intervals")

    linear = IntervalCollection.from_file("graph_reads_on_linear2.intervals").intervals
    #linear = IntervalCollection.create_list_from_file("graph_reads_on_linear.intervals")
    f = open("graph_reads_on_linear.bed", "w")
    path = path.to_indexed_interval()
    linear_reads = []
    for read in linear:
        read.graph = ob_graph
        assert np.all(np.array(read.region_paths) > 0) or np.all(np.array(read.region_paths) < 0)

        dir = "+"
        if read.region_paths[0] < 0:
            dir = "-"
            read = read.get_reverse()

        graph_start = read.start_position
        graph_end = read.end_position

        linear_start = MHC_REGION.start + path.get_offset_at_position(graph_start)
        linear_end = MHC_REGION.start + path.get_offset_at_position(graph_end)

        f.writelines("chr6\t%d\t%d\t.\t0\t%s\n" % (linear_start, linear_end, dir))
    f.close()


