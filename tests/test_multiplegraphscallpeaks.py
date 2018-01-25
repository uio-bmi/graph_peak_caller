from graph_peak_caller.multiplegraphscallpeaks import MultipleGraphsCallpeaks
from offsetbasedgraph import GraphWithReversals as Graph, DirectedInterval, IntervalCollection, Block
import unittest
from graph_peak_caller.snarls import SnarlGraph, SnarlGraphBuilder, SimpleSnarl
from graph_peak_caller.linearsnarls import LinearSnarlMap
from pyvg.sequences import SequenceRetriever


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

        self._create_data()

    def _create_data(self):
        node_offset = 1
        graphs = []
        for chrom_number, chromosome in enumerate(self.chromosomes):

            graph  = Graph({i + node_offset: Block(10) for i in range(0, 3)},
                      {i+node_offset: [i+1+node_offset] for i in range(0, 2)})

            snarls = {
                4+node_offset: SimpleSnarl(node_offset, node_offset+2, 4+node_offset)
            }
            snarlgraph = SnarlGraphBuilder(graph.copy(), snarls,
                                                id_counter=5+node_offset).build_snarl_graphs()
            linear_map = LinearSnarlMap.from_snarl_graph(snarlgraph, graph)
            linear_map_file_name = "test_linear_map_%s.tmp" % chromosome
            linear_map.to_json_files(linear_map_file_name)
            self.linear_maps.append(linear_map_file_name)
            self.sequence_retrievers.append(
                SequenceRetriever({i+node_offset: "A" * 10 for i in range(0, 3)})
            )
            self._create_reads(chrom_number, chromosome, graph)
            node_offset += 3
            graph.to_file(chromosome)

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
                right_sub = peak.get_subinterval(self.fragment_length - self.read_length,
                                                 self.fragment_length)
                right_sub_reverse = right_sub.get_reverse()
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
            sequence_retrievers=self.sequence_retrievers
        )
        self.do_asserts()

    def do_asserts(self):
        for i, chromosome in enumerate(self.chromosomes):
            final_peaks = IntervalCollection.create_list_from_file("multigraphs_" + chromosome + "_max_paths.intervalcollection")
            for peak in self.peaks[i]:
                assert peak in final_peaks


if __name__ == "__main__":
    unittest.main()

