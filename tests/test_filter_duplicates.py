import unittest
from graph_peak_caller.pileup import Pileup
from offsetbasedgraph import Graph, Block, Interval, IntervalCollection
from graph_peak_caller.callpeaks import ExperimentInfo
from graph_peak_caller.sampleandcontrolcreator import SampleAndControlCreator
import numpy as np

class DummyLinearMap:
    def __init__(self):
        pass


class TestFilterDulicates(unittest.TestCase):
    def test_filter_duplicates(self):
        graph = Graph({
                        1: Block(10),
                        2: Block(10),
                        3: Block(10),
                        4: Block(10)
                    },
            {
                1: [2, 4],
                2: [3],
                4: [4]
            })
        graph.to_file("test_graph.tmp")

        intervals = [
            Interval(0, 10, [1, 2, 3]),
            Interval(1, 10, [1, 2, 3]),
            Interval(0, 10, [1, 2, 3])
        ]

        experiment_info = ExperimentInfo(100, 10, 2)

        interval_collection = IntervalCollection(intervals)
        interval_collection.to_file("test_intervals.tmp")

        caller = SampleAndControlCreator("test_graph.tmp", "test_intervals.tmp",
                                     "test_intervals.tmp",
                                     experiment_info=experiment_info,
                                     linear_map="dummy_linear_map")
        filtered_intervals_file_name = caller.filter_duplicates_and_count_intervals(
                        caller.sample_intervals,
                        write_to_file = "filtered_intervals_test.tmp")

        intervals_filtered = []
        for interval in IntervalCollection.from_file(filtered_intervals_file_name, text_file=False):
            intervals_filtered.append(interval)

        self.assertEqual(len(intervals_filtered), len(intervals)-1)
        self.assertEqual(intervals_filtered[0], intervals[0])
        self.assertEqual(intervals_filtered[1], intervals[1])


if __name__ == "__main__":
    unittest.main()