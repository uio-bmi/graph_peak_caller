from graph_peak_caller.multiplegraphscallpeaks import MultipleGraphsCallpeaks
from offsetbasedgraph import GraphWithReversals as Graph, \
    DirectedInterval, IntervalCollection, Block
import unittest
from graph_peak_caller.control.snarls import SnarlGraphBuilder, SimpleSnarl
from graph_peak_caller.control.linearmap import LinearMap
from pyvg.sequences import SequenceRetriever
import logging
import os
logging.basicConfig(level=logging.WARNING)


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

    def _create_data(self):
        node_offset = 1
        for chrom_number, chromosome in enumerate(self.chromosomes):
            graph = Graph(
                {i + node_offset: Block(10) for i in range(0, 3)},
                {i+node_offset: [i+1+node_offset] for i in range(0, 2)})

            linear_map = LinearMap.from_graph(graph)
            linear_map_file_name = "test_linear_map_%s.npz" % chromosome
            linear_map.to_file(linear_map_file_name)
            self.linear_maps.append(linear_map_file_name)
            self.sequence_retrievers.append(
                SequenceRetriever({i+node_offset: "A" * 10
                                   for i in range(0, 3)})
            )
            self._create_reads(chrom_number, chromosome, graph)
            node_offset += 3
            graph.to_file(chromosome + ".obg")

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
        self.sample_reads.append(IntervalCollection(sample_reads))
        self.control_reads.append(IntervalCollection(control_reads))

    def test_run_from_init(self):

        caller = MultipleGraphsCallpeaks(
            self.chromosomes,
            self.chromosomes,
            self.sample_reads,
            self.control_reads,
            self.linear_maps,
            self.fragment_length,
            self.read_length,
            has_control=False,
            sequence_retrievers=None,
            skip_filter_duplicates=True
        )
        caller.run()
        self.do_asserts()

    def test_run_from_init_in_two_steps(self):

        caller = MultipleGraphsCallpeaks(
            self.chromosomes,
            self.chromosomes,
            self.sample_reads,
            self.control_reads,
            self.linear_maps,
            self.fragment_length,
            self.read_length,
            has_control=False,
            sequence_retrievers=None,
            skip_filter_duplicates=True,
            stop_after_p_values=True
        )
        caller.run()

        for i, chromosome in enumerate(self.chromosomes):
            print("Running chrom %s" % chromosome)
            caller = MultipleGraphsCallpeaks(
                self.chromosomes,
                self.chromosomes,
                None,
                None,
                None,
                self.fragment_length,
                self.read_length,
                has_control=False,
                sequence_retrievers=None,
                skip_filter_duplicates=True
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


if __name__ == "__main__":
    unittest.main()

