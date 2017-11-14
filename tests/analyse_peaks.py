from collections import defaultdict
import matplotlib.pyplot as plt
import os
import pyvg
from pyvg.sequences import SequenceRetriever
from pyvg.util import vg_gam_file_to_interval_list
import offsetbasedgraph as obg
from offsetbasedgraph.graphtraverser import GraphTraverserUsingSequence
from graph_peak_caller.peakcollection import PeakCollection
from graph_peak_caller.util import get_linear_paths_in_graph
from offsetbasedgraph import IntervalCollection, DirectedInterval
from graph_peak_caller.subgraphcollection import SubgraphCollection
from graph_peak_caller.peakscores import MaxPathPeakCollection
from .peakscomparer import PeaksComparer, get_peaks_comparer_for_linear_and_graph_peaks


def create_linear_peaks_from_bed(linear_sequence_fasta_file, peaks_bed_file,
                                 obg_graph_file_name, vg_graph_file_name, start_node,
                                 graph_start_offset, graph_end_offset):

    ob_graph = obg.GraphWithReversals.from_file(obg_graph_file_name)
    search_sequence = open(linear_sequence_fasta_file).read()
    sequence_retriever = SequenceRetriever.from_vg_graph(vg_graph_file_name)
    traverser = GraphTraverserUsingSequence(ob_graph, search_sequence, sequence_retriever)
    traverser.search_from_node(start_node)
    linear_path_interval = traverser.get_interval_found()
    IntervalCollection([linear_path_interval]).to_file("linear_path", text_file=True)
    print("Length")
    print(linear_path_interval.length())
    print(linear_path_interval.region_paths[0])
    print(linear_path_interval.start_position)
    print(linear_path_interval.end_position)


    linear_peaks = PeakCollection.create_from_linear_intervals_in_bed_file(
                        obg_graph_file_name,
                        linear_path_interval,
                        peaks_bed_file,
                        graph_start_offset,
                        graph_end_offset)

    linear_peaks.to_file("linear_peaks", text_file=True)


class SubgraphAnalyser(object):

    def __init__(self, graph, subgraphs_file_name):
        self.graph = graph
        self.subgraphs = SubgraphCollection.from_pickle(subgraphs_file_name, graph=graph)



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
        self.subgraphs = SubgraphCollection.from_pickle(subgraphs_file_name, graph=graph)
        self.peaks = PeakCollection.create_list_from_file(peaks_file_name, graph=graph)

    def check_peaks_in_subgraphs(self):
        n_in_subgraphs = 0
        for peak in self.peaks:
            print(peak)
            if self.subgraphs.contains_interval(peak):
                n_in_subgraphs += 1

        print("%d peaks are in subgraphs" % n_in_subgraphs)


class AlignmentsAnalyser(object):
    def __init__(self, vg_graph, vg_gam_file_name, ob_graph, linear_path_interval_file_name):
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
                #if path.contains_in_order_any_direction(read):
                if path.contains(read):
                    hits[path.name] += 1

            if i % 100 == 0:
                print("#%d" % i)
            i += 1

        print("N reads: %d" % i)
        for name, hit in hits.items():
            print("%s: %d " % (name, hit))


#comparer = PeaksComparer("CTCF_peaks.narrowPeak", "real_data_max_paths")
#comparer.compare_q_values_for_similar_peaks()


#comparer.plot_peak_lengths()
#comparer.check_similarity()
#comparer.check_overlap_with_linear_path()
#peaks = comparer.get_peaks_not_on_linear_path()
#comparer.peaks_to_fasta(peaks)




#comparer = SubgraphComparer(ob_graph, "real_data_peaks_as_subgraphs.pickle", "linear_peaks")
#comparer.check_peaks_in_subgraphs()

#analyser = SubgraphAnalyser(ob_graph, "real_data_peaks_as_subgraphs.pickle")
#analyser.print_sizes()

#peaks = comparer.graph_peaks_on_main_path_not_in_linear()
#comparer.peaks_to_fasta(peaks, "alone_linear.peaks")


#analyser = AlignmentsAnalyser(vg_graph, "ENCFF001HNI_filtered_q60.gam", ob_graph, "linear_path")  # sample reads
#analyser = AlignmentsAnalyser(vg_graph, "ENCFF001HNS_filtered_q60.gam", ob_graph, "linear_path")  # Control reads
#analyser.count_alignments_on_linear_paths()
#analyser.count_alignments_on_linear_path()

#compare_linear_and_graph_peaks(ob_graph, "linear_peaks", "real_data_max_paths")
#create_linear_peaks_from_bed("mhc_cleaned2.fa", "../ENCFF155DHA.bed", "cactus-mhc.obg", "cactus-mhc.vg", 225518, 28510119, 33480577)

#graph_peaks = PeakCollection.from_file("real_data_max_paths")

