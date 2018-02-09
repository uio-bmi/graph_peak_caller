from graph_peak_caller.peakcollection import PeakCollection
from graph_peak_caller.nongraphpeaks import NonGraphPeakCollection
from offsetbasedgraph import IntervalCollection
from .analyse_peaks import create_linear_path, LinearRegion
import pyvg


class AnalysisResults:
    def __init__(self):



class PeaksComparerV2(object):

    def __init__(self, graph, vg_json_file_name,
                 linear_peaks_fasta_file_name,
                 graph_peaks_fasta_file_name,
                 linear_peaks_fimo_results_file,
                 graph_peaks_fimo_results_file,
                 region=None):

        self.graph = graph
        self.graph_peaks_fasta_file_name = graph_peaks_fasta_file_name
        self.linear_peaks_fasta_file_name = linear_peaks_fasta_file_name

        self.linear_matching_motif = self._get_peaks_matching_motif(
                linear_peaks_fimo_results_file)
        self.graph_matching_motif = self._get_peaks_matching_motif(
            graph_peaks_fimo_results_file)

        vg_graph = pyvg.vg.Graph.create_from_file(vg_json_file_name)
        self.linear_path = create_linear_path(graph, vg_graph)
        self.peaks1 = PeakCollection.from_fasta_file(self.graph_peaks_fasta_file_name,
                                                     graph)

        nongraphpeaks = NonGraphPeakCollection.from_fasta(self.linear_peaks_fasta_file_name)

        if region is not None:
            nongraphpeaks.filter_peaks_outside_region(region.chrom, region.start, region.end)
        self.peaks2 = PeakCollection.create_from_nongraph_peak_collection(
            graph,
            nongraphpeaks,
            self.linear_path,
            graph_region=region)


        self.results = AnalysisResults()

        self.peaks1_not_in_peaks2 = []
        self.peaks2_not_in_peaks1 = []
        self.run_all_analysis()

    def run_all_analysis(self):
        self.check_similarity()
        self.check_non_matching_for_motif_hits()

    def check_non_matching_for_motif_hits(self):
        print("\n--- Checking peaks that are not in other set for motif hits --- ")
        i = 1
        for not_matching in [self.peaks1_not_in_peaks2, self.peaks2_not_in_peaks1]:
            if i == 1:
                print("Graph peaks")
                matching_motif = self.graph_matching_motif
            else:
                print("Linear peaks")
                matching_motif = self.linear_matching_motif

            n = self.count_matching_motif(matching_motif, not_matching)
            print("N matching motif: %d" % n)
            i += 1

    def count_matching_motif(self, matching_motif_ids, peak_list):
        n = 0
        for peak in peak_list:
            print(peak.unique_id)
            assert peak.unique_id is not None
            if peak.unique_id in matching_motif_ids:
                n += 1
        return n

    def _get_peaks_matching_motif(self, fimo_file_name):
        matches = set()
        f = open(fimo_file_name)
        for line in f:
            l = line.split()
            sequence_id = l[2]
            matches.add(sequence_id)
        f.close()
        return matches

    def _read_peaks(self):
        self.peaks1 = PeakCollection.from_fasta_file(self.graph_peaks_fasta_file_name)

    def plot_peak_lengths(self):
        import matplotlib.pyplot as plt
        i = 1
        for peak_set in (self.peaks1, self.peaks2):
            plt.hist([peak.length() for peak in peak_set],
                     bins=1000, label="Peak set %d" % i)
            i += 1

        plt.legend()
        plt.show()

    def get_peaks_at_same_position(self, allowed_mismatches=0):
        same_pos = []
        for peak in self.peaks1:
            identical = self.peaks2.get_similar_intervals(
                peak, allowed_mismatches)

            assert len(identical) <= 1
            if len(identical) > 0:
                same_pos.append((peak, identical[0]))

        return same_pos

    @classmethod
    def create_from_graph_peaks_and_linear_peaks(
            cls,
            linear_peaks_bed_file_name,
            graph_peaks_file_name,
            ob_graph,
            linear_path,
            graph_region=None):
        linear_path = linear_path.to_indexed_interval()
        linear_peaks = PeakCollection.create_from_linear_intervals_in_bed_file(
            ob_graph, linear_path, linear_peaks_bed_file_name,
            graph_region=graph_region)

        linear_peaks.to_file("macs_peaks.intervals", text_file=True)
        sequence_retriever = None
        comparer = PeaksComparer(ob_graph, sequence_retriever,
                                 None,
                                 "macs_peaks.intervals",
                                 graph_peaks_file_name)
        return comparer

    def check_similarity(self, analyse_first_n_peaks=250):
        i = 1
        for peak_datasets in [(self.peaks1, self.peaks2),
                              (self.peaks2, self.peaks1)]:
            n_identical = 0
            tot_n_similar = 0
            n_similar = 0
            n_tot = 0
            print("\n-- Comparing set %d against set %d ---" % (i, i % 2 + 1))
            peaks1, peaks2 = peak_datasets
            print("Number of peaks in main set: %d" % len(peaks1.intervals))
            not_matching = []
            matching = []
            counter = 0
            for peak in sorted(peaks1, key=lambda x: x.score, reverse=True)[0:analyse_first_n_peaks]:
                assert peak.unique_id is not None
                counter += 1
                if peaks2.contains_interval(peak):
                    n_identical += 1

                similar_intervals = peaks2.get_overlapping_intervals(peak, 50)
                if len(similar_intervals) > 0:
                    n_similar += 1
                    tot_n_similar += len(similar_intervals)
                    matching.append(peak)
                else:
                    not_matching.append(peak)
                    if i == 1:
                        self.peaks1_not_in_peaks2.append(peak)
                    else:
                        self.peaks2_not_in_peaks1.append(peak)

                n_tot += 1

            self.results.peaks1_not_in_peaks2 = len(self.peaks1_not_in_peaks2)
            self.results.peaks2_not_in_peaks1 = len(self.peaks2_not_in_peaks1)

            not_matching = IntervalCollection(not_matching)
            not_matching.to_file("not_matching_set%d.intervals" % i, text_file=True)
            matching = IntervalCollection(matching)
            matching.to_file("matching_set%d.intervals" % i, text_file=True)

            print("Total peaks in main set: %d" % n_tot)
            print("N identical to peak in other set: %d " % n_identical)
            print("N similar to peak in other set: %d " % n_similar)
            print("Total number of similar hits (counting all hits for each peaks): %d " % tot_n_similar)

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
                    n_inside += 1
                if self.linear_path.contains_in_order_any_direction(peak):
                    n_inside_correct_order += 1

            print("n inside: %d / %d" % (n_inside, len(peaks.intervals)))
            print("n inside correct order: %d / %d" % (
                n_inside_correct_order, len(peaks.intervals)))

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


class PeaksComparer(object):

    def __init__(self, graph, sequence_retriever, linear_path_file_name,
                 peaks1_file_name, peaks2_file_name):
        self.graph = graph
        self.sequence_retriever = sequence_retriever
        self.peaks1 = PeakCollection.create_list_from_file(
            peaks1_file_name, graph=graph)
        self.peaks2 = PeakCollection.create_list_from_file(
            peaks2_file_name, graph=graph)
        print("Number of intervals in set 1/2: %d / %d" % (
            len(self.peaks1.intervals), len(self.peaks2.intervals)))
        if linear_path_file_name is not None:
            self.linear_path = IntervalCollection.create_list_from_file(
                linear_path_file_name, self.graph).intervals[0]

    def plot_peak_lengths(self):
        import matplotlib.pyplot as plt
        i = 1
        for peak_set in (self.peaks1, self.peaks2):
            plt.hist([peak.length() for peak in peak_set],
                     bins=1000, label="Peak set %d" % i)
            i += 1

        plt.legend()
        plt.show()

    def get_peaks_at_same_position(self, allowed_mismatches=0):
        same_pos = []
        for peak in self.peaks1:
            identical = self.peaks2.get_similar_intervals(
                peak, allowed_mismatches)

            assert len(identical) <= 1
            if len(identical) > 0:
                same_pos.append((peak, identical[0]))

        return same_pos

    def compare_q_values_for_similar_peaks(self):

        for peak in self.peaks1:
            similar = self.peaks2.get_similar_intervals(
                peak, allowed_mismatches=10)
            if len(similar) > 0:
                print("Found match(es) for %s" % peak)
                for matched_peak in similar:
                    print("   Match agsinst %s with scores %.3f, %.3f" %
                          (matched_peak, peak.score, matched_peak.score))
            else:
                print("No match for peak %s" % peak)

    @classmethod
    def create_from_graph_peaks_and_linear_peaks(
            cls,
            linear_peaks_bed_file_name,
            graph_peaks_file_name,
            ob_graph,
            linear_path,
            graph_region=None):
        linear_path = linear_path.to_indexed_interval()
        linear_peaks = PeakCollection.create_from_linear_intervals_in_bed_file(
            ob_graph, linear_path, linear_peaks_bed_file_name,
            graph_region=graph_region)

        linear_peaks.to_file("macs_peaks.intervals", text_file=True)
        sequence_retriever = None
        comparer = PeaksComparer(ob_graph, sequence_retriever,
                                 None,
                                 "macs_peaks.intervals",
                                 graph_peaks_file_name)
        return comparer

    def check_similarity(self, analyse_first_n_peaks=250):
        i = 1
        for peak_datasets in [(self.peaks1, self.peaks2),
                              (self.peaks2, self.peaks1)]:
            n_identical = 0
            tot_n_similar = 0
            n_similar = 0
            n_tot = 0
            print("\n-- Comparing set %d against set %d ---" % (i, i % 2 + 1))
            peaks1, peaks2 = peak_datasets
            print("Number of peaks in main set: %d" % len(peaks1.intervals))
            not_matching = []
            matching = []
            counter = 0
            for peak in sorted(peaks1, key=lambda x: x.score, reverse=True)[0:analyse_first_n_peaks]:
                counter += 1
                if peaks2.contains_interval(peak):
                    n_identical += 1

                similar_intervals = peaks2.get_overlapping_intervals(peak, 50)
                if len(similar_intervals) > 0:
                    n_similar += 1
                    tot_n_similar += len(similar_intervals)
                    print(counter, peak.score, peak.start_position)
                    matching.append(peak)
                else:
                    not_matching.append(peak)
                    print(peak, "\t", 0)

                n_tot += 1

            not_matching = IntervalCollection(not_matching)
            not_matching.to_file("not_matching_set%d.intervals" % i, text_file=True)
            matching = IntervalCollection(matching)
            matching.to_file("matching_set%d.intervals" % i, text_file=True)

            print("Total peaks in main set: %d" % n_tot)
            print("N identical to peak in other set: %d " % n_identical)
            print("N similar to peak in other set: %d " % n_similar)
            print("Total number of similar hits (counting all hits for each peaks): %d " % tot_n_similar)

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
                    n_inside += 1
                if self.linear_path.contains_in_order_any_direction(peak):
                    n_inside_correct_order += 1

            print("n inside: %d / %d" % (n_inside, len(peaks.intervals)))
            print("n inside correct order: %d / %d" % (
                n_inside_correct_order, len(peaks.intervals)))

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
