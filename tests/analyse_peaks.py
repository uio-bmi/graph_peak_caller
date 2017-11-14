from collections import defaultdict
import matplotlib.pyplot as plt
import pyvg
from pyvg.sequences import SequenceRetriever
from pyvg.util import vg_gam_file_to_interval_list
import offsetbasedgraph as obg
from offsetbasedgraph.graphtraverser import GraphTraverserUsingSequence
from graph_peak_caller.peakcollection import PeakCollection
from graph_peak_caller.util import get_linear_paths_in_graph
from offsetbasedgraph import IntervalCollection, DirectedInterval
from graph_peak_caller.subgraphcollection import SubgraphCollection


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



class PeaksComparer(object):

    def __init__(self, graph, sequence_retriever, linear_path_file_name, peaks1_file_name, peaks2_file_name):
        self.graph = graph
        self.sequence_retriever = sequence_retriever
        self.peaks1 = PeakCollection.create_list_from_file(peaks1_file_name, graph=graph)
        self.peaks2 = PeakCollection.create_list_from_file(peaks2_file_name, graph=graph)
        print("Number of intervals in set 1/2: %d / %d" % (len(self.peaks1.intervals), len(self.peaks2.intervals)))

        self.linear_path = IntervalCollection.create_list_from_file(linear_path_file_name, self.graph).intervals[0]


    def plot_peak_lengths(self):
        import matplotlib.pyplot as plt
        i = 1
        for peak_set in (self.peaks1, self.peaks2):
            plt.hist([peak.length() for peak in peak_set], bins=1000, label="Peak set %d" % i)
            i += 1

        plt.legend()
        plt.show()

    def compare_q_values_for_similar_peaks(self):

        for peak in self.peaks1:
            similar = self.peaks2.get_similar_intervals(peak, allowed_mismatches=10)
            if len(similar) > 0:
                print("Found match(es) for %s" % peak)
                for matched_peak in similar:
                    print("   Match agsinst %s with scores %.3f, %.3f" %
                          (matched_peak, peak.score, matched_peak.score))
            else:
                print("No match for peak %s" % peak)

    def check_similarity(self):
        i = 1
        for peak_datasets in [(self.peaks1, self.peaks2), (self.peaks2, self.peaks1)]:
            n_identical = 0
            tot_n_similar = 0
            n_similar = 0
            n_tot = 0
            print("\n-- Comparing set %d against set %d ---" % (i, i % 2 + 1))
            peaks1, peaks2 = peak_datasets
            print("Number of peaks in main set: %d" % len(peaks1.intervals))
            not_matching = []

            for peak in peaks1:
                if peaks2.contains_interval(peak):
                    n_identical += 1

                similar_intervals = peaks2.get_overlapping_intervals(peak, 50)
                #for similar in similar_intervals:
                #    print("%s is simmilar to %s" % (peak, similar))

                if len(similar_intervals) > 0:
                    #print("%s \n overlaps with \n %s \n\n" % (peak, similar_intervals[0]))
                    n_similar += 1
                    tot_n_similar += len(similar_intervals)
                else:
                    not_matching.append(peak)

                n_tot += 1

            not_matching = IntervalCollection(not_matching)
            not_matching.to_file("not_matching_set%d.intervals" % i)

            print("Total peaks in main set: %d" % n_tot)
            print("N identical to peak in other set: %d " % n_identical)
            print("N similar to peak in other set: %d " % n_similar)
            print("Total number of simmilar hits (counting all hits for each peaks): %d " % tot_n_similar)

            i += 1

    def check_overlap_with_linear_path(self):
        i = 0
        for peaks in [self.peaks1, self.peaks2]:
            print("\n -- Checking peak dataset % against linear path -- " % (i))
            i += 1
            n_inside = 0
            n_inside_correct_order = 0
            for peak in peaks:
                if self.linear_path.contains(peak):
                    #print(peak)
                    n_inside += 1
                if self.linear_path.contains_in_order_any_direction(peak):
                    n_inside_correct_order += 1

                #overlap = linear_interval.overlap(peak)
                #print("Overlap: %d / %d" % (overlap, peak.length()))

            print("n inside: %d / %d" % (n_inside, len(peaks.intervals)))
            print("n inside correct order: %d / %d" % (n_inside_correct_order, len(peaks.intervals)))

    def get_peaks_on_linear_path(self, peaks):
        out = []
        for peak in peaks:
            if self.linear_path.contains_in_correct_order(peak):
                out.append(peak)
        return out

    def get_peaks_not_on_linear_path(self):
        out = []
        for peak in self.peaks2:
            if not self.linear_path.contains_in_correct_order(peak):
                out.append(peak)
        return out

    def peaks_to_fasta(self, peaks, file_name="peaks.fasta"):
        f = open(file_name, "w")
        i = 1
        for p in peaks:
            sequence = self.sequence_retriever.get_interval_sequence(p)
            f.writelines([">peak%d\n%s\n" % (i, sequence)])
            i += 1
        f.close()

    def graph_peaks_on_main_path_not_in_linear(self):
        on_linear = self.get_peaks_on_linear_path(self.peaks2)
        print("On linear path: %d" % len(on_linear))
        on_linear = PeakCollection(on_linear)
        return on_linear.get_peaks_not_in_other_collection(self.peaks1, 20)


class AlignmentsAnalyser(object):
    def __init__(self, vg_graph, vg_gam_file_name, ob_graph, linear_path_interval_file_name):
        self.graph = ob_graph
        self.vg_graph = vg_graph
        print("Reading reads")
        self.reads = vg_gam_file_to_interval_list(vg_graph, vg_gam_file_name, ob_graph, max_intervals=10000)
        print("Number of reads: %d" % len(self.reads))

        self.linear_path = IntervalCollection.create_list_from_file(linear_path_interval_file_name, self.graph).intervals[0]

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



def get_peaks_comparer_for_linear_and_graph_peaks(
                linear_peaks_bed_file_name,
                graph_peaks_file_name):
    import os
    sequence_retriever = None  # SequenceRetriever.from_vg_graph("haplo1kg50-mhc.vg")
    ob_graph = obg.GraphWithReversals.from_file("haplo1kg50-mhc.obg")
    if not os.isfile("linear_paths_haplo1kg-50.intervals"):
        vg_graph = pyvg.Graph.create_from_file("haplo1kg50-mhc.json")
        linear_paths = get_linear_paths_in_graph(
                        ob_graph,
                        vg_graph,
                        write_to_file_name="linear_paths_haplo1kg-50.intervals")

    linear_path = IntervalCollection.create_list_from_file(
                        "linear_paths_haplo1kg-50.intervals",
                        ob_graph).intervals[0]
    linear_path = linear_path.to_indexed_interval()
    linear_peaks = PeakCollection.create_from_linear_intervals_in_bed_file(ob_graph,
                                                                           linear_path,
                                                                           linear_peaks_bed_file_name,
                                                                           28510119,
                                                                           33480577)
    linear_peaks.to_file("mac_peaks.intervals", text_file=True)
    comparer = PeaksComparer(ob_graph, sequence_retriever,
                             "linear_paths_haplo1kg-50.intervals",
                             "mac_peaks.intervals",
                             graph_peaks_file_name)


comparer = get_peaks_comparer_for_linear_and_graph_peaks(
            "CTCF_peaks.narrowPeak", "real_data_max_paths")

comparer.compare_q_values_for_similar_peaks()


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

