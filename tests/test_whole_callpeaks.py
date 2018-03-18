import unittest
import logging
from offsetbasedgraph import GraphWithReversals, Block, \
    DirectedInterval, IntervalCollection
from graph_peak_caller import ExperimentInfo, CallPeaks, Configuration
from graph_peak_caller.control.linearmap import LinearMap
from graph_peak_caller.reporter import Reporter
from graph_peak_caller.intervals import Intervals
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s, %(levelname)s: %(message)s")


class TestWholeCallPeaks(unittest.TestCase):
    def _create_reads(self):
        self.sample_reads = []
        for peak in self.peaks:
            for i in range(0, 10):
                left_sub = peak.get_subinterval(0, self.read_length)
                self.sample_reads.append(left_sub)
                right_sub = peak.get_subinterval(
                    self.fragment_length - self.read_length,
                    self.fragment_length)
                right_sub_reverse = right_sub.get_reverse()
                self.sample_reads.append(right_sub_reverse)

    def assert_final_peaks_equals_input_peaks(self):
        final_peaks = IntervalCollection.create_list_from_file(
            "test_max_paths.intervalcollection")
        for peak in self.peaks:
            self.assertTrue(peak in final_peaks.intervals,
                            "Peak %s not in final peaks. Final peaks: \n%s"
                            % (peak, final_peaks.intervals))
        self.assertEqual(len(self.peaks), len(final_peaks.intervals))

    def setUp(self):
        self.set_graph()

    def run_callpeaks(self):
        self._create_reads()
        control_reads = self.sample_reads.copy()

        self.graph_size = sum(block.length()
                              for block in self.graph.blocks.values())
        config = Configuration()
        config.fragment_length = self.fragment_length
        config.read_length = self.read_length
        config.linear_map_name = "test_linear_map.npz"
        caller = CallPeaks(self.graph, config, Reporter("test_"))
        caller.run(Intervals(self.sample_reads), Intervals(control_reads))

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
        LinearMap.from_graph(self.graph).to_file("test_linear_map.npz")

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
        blocks = {i: Block(3) for i in range(1, 11)}
        blocks[11] = Block(1000)
        self.graph = GraphWithReversals(
            blocks,
            {
                1: [2, 3],
                2: [7, 8],
                3: [4, 5],
                4: [6],
                5: [6],
                6: [10],
                7: [9],
                8: [9],
                9: [10],
                10: [11]
             })

        LinearMap.from_graph(self.graph).to_file("test_linear_map.npz")

    def test_single_peak(self):
        self.peaks = [
            DirectedInterval(1, 1, [1, 2, 7], self.graph),
        ]
        self.do_asserts()

    def test_multiple_peaks(self):
        self.peaks = [
            DirectedInterval(2, 2, [1, 3, 4], self.graph),
            DirectedInterval(2, 2, [6, 10, 11], self.graph),
        ]
        self.do_asserts()

    def test_multiple_peaks2(self):

        self.peaks = [
            DirectedInterval(2, 2, [1, 3, 4], self.graph),
            DirectedInterval(2, 2, [6, 10, 11], self.graph),
        ]
        self.do_asserts()

    def test_multiple_peaks3(self):

        self.peaks = [
            DirectedInterval(2, 2, [1, 3, 4], self.graph),
            DirectedInterval(0, 3, [7, 9], self.graph)
        ]
        self.do_asserts()

if __name__ == "__main__":
    unittest.main()
