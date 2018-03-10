from graph_peak_caller.peakcollection import PeakCollection
from graph_peak_caller.nongraphpeaks import NonGraphPeakCollection
from offsetbasedgraph import IntervalCollection
from offsetbasedgraph import IndexedInterval, NumpyIndexedInterval
from graph_peak_caller.util import create_linear_path
import pickle
import pyvg
import logging


class AnalysisResults:
    def __init__(self):
        self.tot_peaks1 = 0
        self.tot_peaks2 = 0
        self.peaks1_in_peaks2 = 0
        self.peaks2_in_peaks1 = 0
        self.peaks2_not_in_peaks1 = 0
        self.peaks1_not_in_peaks2 = 0
        self.peaks2_not_in_peaks1_matching_motif = 0
        self.peaks1_not_in_peaks2_matching_motif = 0
        self.peaks1_in_peaks2_matching_motif = 0
        self.peaks2_in_peaks1_matching_motif = 0
        self.motif_ambiguous = 0
        self.motif_not_ambiguous = 0
        self.not_motif_ambiguous = 0
        self.not_motif_not_ambiguous = 0

    def __repr__(self):
        out = ""
        out += "Number of peaks1: %d \n" % self.tot_peaks1
        out += "Number of peaks2: %d \n" % self.tot_peaks2
        out += "Number of peaks1 also found in peaks2: %d (%d with motif hit) \n" % \
               (self.peaks1_in_peaks2, self.peaks1_in_peaks2_matching_motif)
        out += "Number of peaks2 also found in peaks1: %d (%d with motif hit) \n" % \
               (self.peaks2_in_peaks1, self.peaks2_in_peaks1_matching_motif)

        out += "\n"
        out += "Number of peaks1 NOT found in peaks 2: %d (%d with motif hit) \n" % \
               (self.peaks1_not_in_peaks2, self.peaks1_not_in_peaks2_matching_motif)
        out += "Number of peaks2 NOT found in peaks 1: %d (%d with motif hit) \n" % \
               (self.peaks2_not_in_peaks1, self.peaks2_not_in_peaks1_matching_motif)
        out += "\n \n"
        return out

    def __str__(self):
        return self.__repr__()

    def __add__(self, other):
        self.tot_peaks1 += other.tot_peaks1
        self.tot_peaks2 += other.tot_peaks2
        self.peaks1_in_peaks2 += other.peaks1_in_peaks2
        self.peaks2_in_peaks1 += other.peaks2_in_peaks1
        self.peaks2_not_in_peaks1 += other.peaks2_not_in_peaks1
        self.peaks1_not_in_peaks2 += other.peaks1_not_in_peaks2
        self.peaks2_not_in_peaks1_matching_motif += other.peaks2_not_in_peaks1_matching_motif
        self.peaks1_not_in_peaks2_matching_motif += other.peaks1_not_in_peaks2_matching_motif
        self.peaks2_in_peaks1_matching_motif += other.peaks2_in_peaks1_matching_motif
        self.peaks1_in_peaks2_matching_motif += other.peaks1_in_peaks2_matching_motif
        self.motif_ambiguous += other.motif_ambiguous
        self.motif_not_ambiguous += other.motif_not_ambiguous
        self.not_motif_ambiguous += other.not_motif_ambiguous
        self.not_motif_not_ambiguous += other.not_motif_not_ambiguous
        return self

    def to_file(self, file_name):
        with open(file_name, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def from_file(cls, file_name):
        f = open(file_name, "rb")
        results = pickle.load(f)
        f.close()
        return results

# $ graph_peak_caller analyse_peaks graph.nobg haplo1kg50-mhc.json  macs_sequences_mhc.fasta test_sequences.fasta  fimo_macs_sequences/fimo.txt fimo_test_sequences/fimo.txt chr6 28510119 33480577
# $ graph_peak_caller analyse_peaks ~/dev/graph_peak_caller/tests/whole_genome/22.nobg ~/dev/graph_peak_caller/tests/whole_genome/22.json macs_sequences_chr22.fasta 22_sequences.fasta fimo_macs_chr22/fimo.txt fimo_graph_chr22/fimo.txt 22 22 0 0


class PeaksComparerV2(object):
    def __init__(self, graph,
                 linear_peaks_fasta_file_name,
                 graph_peaks_fasta_file_name,
                 linear_peaks_fimo_results_file,
                 graph_peaks_fimo_results_file,
                 linear_path,
                 region=None):

        assert isinstance(linear_path, NumpyIndexedInterval), \
            "Linear path should be numpy indexed interval for fast comparison"

        print(linear_path)

        self.graph = graph
        self.graph_peaks_fasta_file_name = graph_peaks_fasta_file_name
        self.linear_peaks_fasta_file_name = linear_peaks_fasta_file_name

        self.linear_matching_motif = self._get_peaks_matching_motif(
                linear_peaks_fimo_results_file)
        self.graph_matching_motif = self._get_peaks_matching_motif(
            graph_peaks_fimo_results_file)

        self.linear_path = linear_path
        self.peaks1 = PeakCollection.from_fasta_file(self.graph_peaks_fasta_file_name,
                                                     graph)

        self.peaks1.create_node_index()

        nongraphpeaks = NonGraphPeakCollection.from_fasta(self.linear_peaks_fasta_file_name)

        if region is not None:
            nongraphpeaks.filter_peaks_outside_region(
                region.chromosome, region.start, region.end)
        else:
            logging.info(
                "Graph region is None. Whole graph will be used and no filtering on peaks will be done.")

        self.peaks2 = PeakCollection.create_from_nongraph_peak_collection(
            graph,
            nongraphpeaks,
            self.linear_path,
            graph_region=region)

        print("%d linear peaks converted to graph" % len(self.peaks2.intervals))

        self.peaks2.create_node_index()

        self.results = AnalysisResults()

        self.peaks1_in_peaks2 = []
        self.peaks2_in_peaks1 = []
        self.peaks1_not_in_peaks2 = []
        self.peaks2_not_in_peaks1 = []
        self.run_all_analysis()

    def run_all_analysis(self):
        self.check_similarity()
        self.check_non_matching_for_motif_hits()
        self.check_matching_for_motif_hits()

        print(self.results)

    def check_matching_for_motif_hits(self):
        print("\n--- Checking peaks matching for motif hits --- ")
        self._get_pair_ambig_status()
        i = 1
        for matching in [self.peaks1_in_peaks2, self.peaks2_in_peaks1]:
            if i == 1:
                print("Checking graph peaks, in total %d in other set" % len(matching))
                matching_motif = self.graph_matching_motif
            else:
                print("Checking linear peaks, in total %d in other set" % len(matching))
                matching_motif = self.linear_matching_motif

            a, b, c, d = self.count_matching_motif(
                matching_motif, matching, with_info=True)
            n = a + b
            print("Found n peaks that matches motif: %d" % n)

            if i == 1:
                self.results.motif_ambiguous = a
                self.results.motif_not_ambiguous = b
                self.results.not_motif_ambiguous = c
                self.results.not_motif_not_ambiguous = d
                self.results.peaks1_in_peaks2_matching_motif = n
            else:
                print("MA: ", a, "MNA:", b, "NMA:", c, "NMNA:", d)
                self.results.peaks2_in_peaks1_matching_motif = n
            i += 1

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

            if i == 1:
                self.results.peaks1_not_in_peaks2_matching_motif = n
            else:
                self.results.peaks2_not_in_peaks1_matching_motif = n

            i += 1

    def _get_pair_ambig_status(self):
        peaks = self.peaks1_in_peaks2
        matches = [self.peaks2.which_approx_contains_part_of_interval(peak)
                   for peak in peaks]
        my_motif_matches = [peak.unique_id in self.graph_matching_motif
                            for peak in peaks]
        other_motif_matches = [peak.unique_id in self.linear_matching_motif
                               for peak in peaks]
        is_ambiguous = [peak.info[0] for peak in peaks]
        n = 0
        found = 0
        for my_match, other_match, ambig in zip(my_motif_matches, other_motif_matches, is_ambiguous):
            if other_match and not my_match:
                n += 1
                if ambig:
                    found += 1
        print("Ambigous of non-found: ", found, n)

    def count_matching_motif(self, matching_motif_ids, peak_list, with_info=False):
        a, b, c, d = (0, 0, 0, 0)
        for peak in peak_list:
            assert peak.unique_id is not None
            ambig = peak.info[0]
            motif = peak.unique_id in matching_motif_ids
            if motif and ambig:
                a += 1
            elif motif and not ambig:
                b += 1
            elif (not motif) and ambig:
                c += 1
            else:
                d += 1
        if with_info:
            return a, b, c, d
        return a+b

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

    def check_similarity(self, analyse_first_n_peaks=10000000):
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
            if i == 1:
                self.results.tot_peaks1 = len(peaks1.intervals)
            else:
                self.results.tot_peaks2 = len(peaks1.intervals)

            not_matching = []
            matching = []
            counter = 0
            for peak in sorted(peaks1, key=lambda x: x.score, reverse=True)[0:analyse_first_n_peaks]:
                assert peak.unique_id is not None
                counter += 1
                #if peaks2.contains_interval(peak):
                #    n_identical += 1

                if counter % 500 == 0:
                    logging.info("Checked %d peaks" % counter)


                #similar_intervals = peaks2.get_overlapping_intervals(peak, 50)
                #print("Peak %s is overlapping with %d other peaks" % (peak.unique_id, len(similar_intervals)))
                #if len(similar_intervals) > 0:
                if peaks2.approx_contains_part_of_interval(peak):
                    n_similar += 1
                    if i == 1:
                        self.peaks1_in_peaks2.append(peak)
                    else:
                        self.peaks2_in_peaks1.append(peak)

                    matching.append(peak)
                else:
                    not_matching.append(peak)
                    if i == 1:
                        self.peaks1_not_in_peaks2.append(peak)
                    else:
                        self.peaks2_not_in_peaks1.append(peak)

                n_tot += 1
            self.results.peaks1_in_peaks2 = len(self.peaks1_in_peaks2)
            self.results.peaks2_in_peaks1 = len(self.peaks2_in_peaks1)

            self.results.peaks1_not_in_peaks2 = len(self.peaks1_not_in_peaks2)
            self.results.peaks2_not_in_peaks1 = len(self.peaks2_not_in_peaks1)

            not_matching = IntervalCollection(not_matching)
            not_matching.to_file("not_matching_set%d.intervals" % i, text_file=True)
            logging.info("Wrote peaks not matching to file not_matching_set%d.intervals" % i)
            matching = IntervalCollection(matching)
            matching.to_file("matching_set%d.intervals" % i, text_file=True)

            logging.info("Total peaks in main set: %d" % n_tot)
            logging.info("N similar to peak in other set: %d " % n_similar)
            logging.info("N not matching other set: %d " % len(not_matching.intervals))

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

            not_matching = PeakCollection(not_matching)
            not_matching.to_file("not_matching_set%d.intervals" % i, text_file=True)
            matching = PeakCollection(matching)
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
