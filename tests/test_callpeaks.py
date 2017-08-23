import unittest
from graph_peak_caller.pileup import Pileup
from offsetbasedgraph import Graph, Block, Interval, IntervalCollection
from graph_peak_caller import bdgcmp, callpeaks
import numpy as np



class TestCallpeaks(unittest.TestCase):
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
        graph.to_file("test_graph")

        intervals = [
            Interval(0, 10, [1, 2, 3]),
            Interval(1, 10, [1, 2, 3]),
            Interval(0, 10, [1, 2, 3])
        ]
        interval_collection = IntervalCollection(intervals)
        interval_collection.to_file("test_intervals")

        caller = callpeaks.CallPeaks("test_graph", "test_intervals")
        filtered_intervals_file_name = caller.filter_duplicates(caller.sample_file_name)

        intervals_filtered = []
        for interval in IntervalCollection.create_generator_from_file(filtered_intervals_file_name):
            intervals_filtered.append(interval)

        self.assertEqual(len(intervals_filtered), len(intervals)-1)
        self.assertEqual(intervals_filtered[0], intervals[0])
        self.assertEqual(intervals_filtered[1], intervals[1])


if __name__ == "__main__":
    unittest.main()