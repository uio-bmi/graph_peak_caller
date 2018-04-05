from graph_peak_caller.multiplegraphscallpeaks import MultipleGraphsCallpeaks
from graph_peak_caller.intervals import Intervals
from graph_peak_caller import Configuration
from graph_peak_caller.reporter import Reporter
from offsetbasedgraph import GraphWithReversals as Graph, \
    DirectedInterval, IntervalCollection, Block, SequenceGraph, Interval
import unittest
from graph_peak_caller.control.linearmap import LinearMap
from pyvg.sequences import SequenceRetriever
import logging
from graph_peak_caller.logging_config import set_logging_config
#set_logging_config(1)
import os
from graph_peak_caller.command_line_interface import run_argument_parser



class TestMultipleGraphsCallPeaks(unittest.TestCase):
    def setUp(self):
        self.chromosomes = ["1", "2", "3", "X", "Y"]
        self.fragment_length = 5
        self.read_length = 2
        self.sample_reads = []
        self.control_reads = []
        self.linear_maps = []
        self.sequence_retrievers = []
        self.peaks = []

        for chrom in self.chromosomes:
            # Delete old files if existing
            if os.path.isfile("multigraphs_%s_pvalues_indexes.npy" % chrom):
                os.remove("multigraphs_%s_pvalues_indexes.npy" % chrom)
                os.remove("multigraphs_%s_pvalues_values.npy" % chrom)

            # Delete old files if existing
            if os.path.isfile("multigraphs_%s_max_paths.intervalcollection" % chrom):
                os.remove("multigraphs_%s_max_paths.intervalcollection" % chrom)

        self._create_data()
        self.config = Configuration()
        self.config.fragment_length = self.fragment_length
        self.config.read_length = self.read_length
        self.config.has_control = False
        self.config.min_background = 0.33
        self.reporter = Reporter("multigraphs_")

    def _create_data(self):
        node_offset = 1
        for chrom_number, chromosome in enumerate(self.chromosomes):
            graph = Graph(
                {i + node_offset: Block(10) for i in range(0, 3)},
                {i+node_offset: [i+1+node_offset] for i in range(0, 2)})

            linear_map = LinearMap.from_graph(graph)
            linear_map_file_name = "linear_map_%s.npz" % chromosome
            linear_map.to_file(linear_map_file_name)
            self.linear_maps.append(linear_map_file_name)
            self.sequence_retrievers.append(
                SequenceRetriever({i+node_offset: "A" * 10
                                   for i in range(0, 3)})
            )
            self._create_reads(chrom_number, chromosome, graph)
            node_offset += 3
            graph.convert_to_numpy_backend()
            SequenceGraph.create_empty_from_ob_graph(graph).to_file(chromosome + ".nobg.sequences")
            graph.to_file(chromosome + ".nobg")

    def _create_reads(self, chrom_number, chrom, graph):
        i = chrom_number
        sample_reads = []
        control_reads = []
        peaks = [DirectedInterval(7, 2, [1 + 3*i, 2 + 3*i], graph)]
        self.peaks.append(peaks)
        for peak in peaks:
            for i in range(0, 10):
                left_sub = peak.get_subinterval(0, self.read_length)
                sample_reads.append(left_sub)
                control_reads.append(left_sub)
                right_sub = peak.get_subinterval(
                    self.fragment_length - self.read_length,
                    self.fragment_length)
                right_sub_reverse = right_sub.get_reverse()
                sample_reads.append(right_sub_reverse)
                control_reads.append(right_sub_reverse)
        self.sample_reads.append(Intervals(sample_reads))
        self.control_reads.append(Intervals(control_reads))

    def test_run_from_init(self):
        caller = MultipleGraphsCallpeaks(
            self.chromosomes,
            [chrom + ".nobg" for chrom in self.chromosomes],
            self.sample_reads,
            self.control_reads,
            self.linear_maps,
            self.config,
            self.reporter
        )
        caller.run()
        self.do_asserts()

    def test_run_from_init_in_two_steps(self):

        caller = MultipleGraphsCallpeaks(
            self.chromosomes,
            [chrom + ".nobg" for chrom in self.chromosomes],
            self.sample_reads,
            self.control_reads,
            self.linear_maps,
            self.config,
            self.reporter,
            stop_after_p_values=True
        )
        caller.run()

        for i, chromosome in enumerate(self.chromosomes):
            caller = MultipleGraphsCallpeaks(
                self.chromosomes,
                 [chrom + ".nobg" for chrom in self.chromosomes],
                None,
                None,
                None,
                self.config,
                self.reporter
            )
            caller.create_joined_q_value_mapping()
            caller.run_from_p_values(only_chromosome=chromosome)
        self.do_asserts()

    def do_asserts(self):
        for i, chromosome in enumerate(self.chromosomes):
            final_peaks = IntervalCollection.create_list_from_file(
                "multigraphs_" + chromosome + "_max_paths.intervalcollection")
            for peak in self.peaks[i]:
                assert peak in final_peaks


class TestMultipleGraphsCallPeaksCommandLine(TestMultipleGraphsCallPeaks):
    # Same test, but using commmand line interface

    def _create_reads(self, *args):
        super(TestMultipleGraphsCallPeaksCommandLine, self)._create_reads(*args)
        for intervals, chrom in zip(self.sample_reads, self.chromosomes):
            IntervalCollection(intervals._intervals).to_file("test_sample_" + chrom + ".intervalcollection", text_file=True)

    def test_typical_run(self):

        print(" ========= Running start ====")
        run_argument_parser(["callpeaks",
                             "-g", "*.nobg",
                             "-s", "test_sample_*.intervalcollection",
                             "-f", "%s" % self.fragment_length,
                             "-r", "%s" % self.read_length,
                             "-u", "100",
                             "-G", "150",
                             "-n", "multigraphs_",
                             "-p", "True",
                             "-D", "True"])


        for i, chromosome in enumerate(self.chromosomes):
            run_argument_parser(["callpeaks_whole_genome_from_p_values", chromosome,
                                 "-d", "./",
                                 "-f", "%s" % self.fragment_length,
                                 "-r", "%s" % self.read_length,
                                 "-n", "multigraphs_"])
        self.do_asserts()

    def test_count_unique_reads(self):
        reads = [
            IntervalCollection([
                Interval(4, 10, [1, 2, 3]),
                Interval(4, 5, [1]),
                Interval(5, 5, [1]),
                Interval(6, 2, [-3, -2, -1])
            ])
        ]
        unique = MultipleGraphsCallpeaks.count_number_of_unique_reads(reads)
        self.assertEqual(unique, 3)




if __name__ == "__main__":
    unittest.main()

