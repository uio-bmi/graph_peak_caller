from graph_peak_caller.peakcollection import PeakCollection
from graph_peak_caller.util import bed_intervals_to_graph
import offsetbasedgraph as obg
from offsetbasedgraph.graphtraverser import GraphTraverserUsingSequence
from pyvg.sequences import SequenceRetriever
from pybedtools import BedTool
from pyvg.util import vg_gam_file_to_interval_collection, vg_gam_file_to_interval_list
from graph_peak_caller.util import get_linear_paths_in_graph
from offsetbasedgraph import IntervalCollection
import pyvg
from collections import defaultdict

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


class PeaksComparer(object):

    def __init__(self, graph, sequence_retriever, linear_path_file_name, peaks1_file_name, peaks2_file_name):
        self.graph = graph
        self.sequence_retriever = sequence_retriever
        self.peaks1 = PeakCollection.create_list_from_file(peaks1_file_name, graph=graph)
        self.peaks2 = PeakCollection.create_list_from_file(peaks2_file_name, graph=graph)
        self.linear_path = IntervalCollection.create_list_from_file(linear_path_file_name, self.graph).intervals[0]

    def check_similarity(self):
        n_identical = 0
        n_similar = 0
        n_tot = 0
        for peak in self.peaks1:
            #print("Linear peak: %s" % peak)
            if self.peaks2.contains_interval(peak):
                print("  IDENTICAL MATCH")
                n_identical += 1

            similar_intervals = self.peaks2.get_similar_intervals(peak, 50)
            #for similar in similar_intervals:
            #    print("%s is simmilar to %s" % (peak, similar))

            if len(similar_intervals) > 0:
                n_similar += 1

            n_tot += 1
        print("Total linear peaks: %d" % n_tot)
        print("N identical: %d " % n_identical)
        print("N similar: %d " % n_similar)

    def check_overlap_with_linear_path(self):

        for peaks in [self.peaks1, self.peaks2]:
            print("Peak dataset ========")
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



"""
sequence_retriever = SequenceRetriever.from_vg_graph("cactus-mhc.vg")
ob_graph = obg.GraphWithReversals.from_file("cactus-mhc.obg")
vg_graph = vg_graph = pyvg.Graph.create_from_file("cactus-mhc.json")
"""

sequence_retriever = SequenceRetriever.from_vg_graph("haplo1kg50-mhc.vg")
ob_graph = obg.GraphWithReversals.from_file("haplo1kg50-mhc.obg")
vg_graph = vg_graph = pyvg.Graph.create_from_file("haplo1kg50-mhc.json")

comparer = PeaksComparer(ob_graph, sequence_retriever, "linear_path", "linear_peaks", "real_data_max_paths")
#comparer.check_similarity()
comparer.check_overlap_with_linear_path()
#peaks = comparer.get_peaks_not_on_linear_path()
#comparer.peaks_to_fasta(peaks)

#peaks = comparer.graph_peaks_on_main_path_not_in_linear()
#comparer.peaks_to_fasta(peaks, "alone_linear.peaks")


#analyser = AlignmentsAnalyser(vg_graph, "ENCFF001HNI_filtered_q60.gam", ob_graph, "linear_path")  # sample reads
#analyser = AlignmentsAnalyser(vg_graph, "ENCFF001HNS_filtered_q60.gam", ob_graph, "linear_path")  # Control reads
#analyser.count_alignments_on_linear_paths()
#analyser.count_alignments_on_linear_path()

#compare_linear_and_graph_peaks(ob_graph, "linear_peaks", "real_data_max_paths")
#create_linear_peaks_from_bed("mhc_cleaned2.fa", "../ENCFF155DHA.bed", "cactus-mhc.obg", "cactus-mhc.vg", 225518, 28510119, 33480577)

#linear_peaks = PeakCollection.create_from_linear_intervals_in_bed_file("..ENCFF155DHA.bed")
#graph_peaks = PeakCollection.from_file("real_data_max_paths")

