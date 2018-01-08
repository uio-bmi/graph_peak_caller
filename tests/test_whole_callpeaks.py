import unittest
from offsetbasedgraph import GraphWithReversals, Block, \
        Interval, DirectedInterval, IntervalCollection
from graph_peak_caller.callpeaks import ExperimentInfo, CallPeaks
from graph_peak_caller.snarls import SnarlGraph, SnarlGraphBuilder, SimpleSnarl
from graph_peak_caller.linearsnarls import LinearSnarlMap
from graph_peak_caller.sparsepileup import SparsePileup
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s, %(levelname)s: %(message)s")

class TestWholeCallPeaks(unittest.TestCase):
    def get_caller(self, sample_intervals, control_intervals, has_control=False):
        self.graph_size = sum(block.length() for block in self.graph.blocks.values())
        experiment_info = ExperimentInfo(self.graph_size,
                                         self.fragment_length, self.read_length)
        return CallPeaks(self.graph,
                         IntervalCollection(sample_intervals),
                         IntervalCollection(control_intervals),
                         experiment_info,
                         out_file_base_name="test_", has_control=has_control,
                         linear_map="test_linear_map.tmp",
                         skip_filter_duplicates=True)

    def _create_reads(self):
        self.sample_reads = []
        for peak in self.peaks:
            for i in range(0, 10):
                left_sub = peak.get_subinterval(0, self.read_length)
                self.sample_reads.append(left_sub)
                right_sub = peak.get_subinterval(self.fragment_length - self.read_length,
                                                 self.fragment_length)
                right_sub_reverse = right_sub.get_reverse()
                self.sample_reads.append(right_sub_reverse)

    def assert_final_peaks_equals_input_peaks(self):
        final_peaks = IntervalCollection.create_list_from_file("test_max_paths.intervalcollection")
        for peak in self.peaks:
            self.assertTrue(peak in final_peaks.intervals, "Peak %s not in final peaks. Final peaks: \n%s" % (peak, final_peaks.intervals))
        self.assertEqual(len(self.peaks), len(final_peaks.intervals))

    def setUp(self):
        self.set_graph()

    def run_callpeaks(self):
        self._create_reads()
        #for read in self.sample_reads:
        #    print(read)

        control_reads = self.sample_reads.copy()

        self.caller = self.get_caller(self.sample_reads,
                                      control_reads,
                                      has_control=False)
        self.caller.run()

    def do_asserts(self):
        for peak in self.peaks:
            assert peak.length() == self.fragment_length

        self.run_callpeaks()
        self.assert_final_peaks_equals_input_peaks()


class TestWholeCallpeaksSplitGraph(TestWholeCallPeaks):

    def set_graph(self):
        self.fragment_length = 5
        self.read_length = 1

        self.graph = GraphWithReversals(
            {i: Block(15) for i in range(1, 5)},
            {
                1: [2, 3],
                2: [4],
                3: [4]
            }
        )
        self.snarls = {
            5: SimpleSnarl(1, 4, 5)
        }
        self.snarlgraph = SnarlGraphBuilder(self.graph.copy(), self.snarls,
                                            id_counter=6).build_snarl_graphs()
        LinearSnarlMap.from_snarl_graph(self.snarlgraph, self.graph).to_json_files("test_linear_map.tmp")

    def test_single_peak(self):
        self.peaks = [
            DirectedInterval(3, 8, [2], self.graph),
        ]
        self.do_asserts()

    def test_peak_crossing_blocks(self):
        self.peaks = [
            DirectedInterval(13, 3, [1, 2], self.graph),
            DirectedInterval(13, 3, [2, 4], self.graph),
        ]
        self.do_asserts()

    def test_multiple_peaks(self):
        self.peaks = [
            DirectedInterval(3, 8, [1], self.graph),
            DirectedInterval(7, 12, [4], self.graph),
            DirectedInterval(14, 4, [3, 4], self.graph),
        ]
        self.do_asserts()


class TestWholeCallPeaksHierarchical(TestWholeCallPeaks):
    def set_graph(self):
        self.fragment_length = 6
        self.read_length = 2

        self.graph = GraphWithReversals(
            {i: Block(3) for i in range(1, 13)},
            {
                11: [1],
                1: [2, 3],
                2: [7, 8],
                3: [4, 5],
                4: [6],
                5: [6],
                6: [10],
                7: [9],
                8: [9],
                9: [10],
                10: [12]
             })

        self.graph.blocks[11] = Block(500)
        self.graph.blocks[12] = Block(500)

        self.snarlgraph = SnarlGraph(
            {
                11: Block(500),
                12: Block(500),
                1: Block(3),
                10: Block(3),
                20: SnarlGraph(
                    {
                        3: Block(3),
                        21: SnarlGraph(
                            {
                                4: Block(3),
                                5: Block(3)
                            },
                            {
                                3: [4, 5],
                                4: [6],
                                5: [6]
                            },
                            start_node=3,
                            end_node=6
                        ),
                        22: SnarlGraph(
                            {
                                7: Block(3),
                                8: Block(3)
                            },
                            {
                                2: [7, 8],
                                7: [9],
                                8: [9]
                            },
                            start_node=2,
                            end_node=9
                        ),
                        2: Block(3),
                        6: Block(3),
                        9: Block(3),
                    },
                    {
                        3: [21],
                        2: [22],
                        21: [6],
                        22: [9],
                        1: [2, 3],
                        6: [10],
                        9: [10]
                    },
                    start_node=1,
                    end_node=10
                )
            },
            {
                11: [1],
                1: [20],
                20: [10],
                10: [12],
                13: [11],  # Dummy
                12: [14],   # Dummy
            },
            start_node=13,
            end_node=14
        )

        LinearSnarlMap.from_snarl_graph(self.snarlgraph, self.graph).to_json_files("test_linear_map.tmp")

    def test_single_peak(self):
        
        self.peaks = [
            DirectedInterval(1, 1, [1, 2, 7], self.graph),
        ]
        self.do_asserts()

    def test_multiple_peaks(self):
        
        self.peaks = [
            DirectedInterval(2, 2, [1, 3, 4], self.graph),
            DirectedInterval(2, 2, [6, 10, 12], self.graph),
        ]
        self.do_asserts()

if __name__ == "__main__":
    unittest.main()
