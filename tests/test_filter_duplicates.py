import unittest
from offsetbasedgraph import Interval, IntervalCollection
from graph_peak_caller.sampleandcontrolcreator import UniqueIntervals


class DummyLinearMap:
    def __init__(self):
        pass


class TestFilterDulicates(unittest.TestCase):
    def test_filter_duplicates(self):
        intervals = [
            Interval(0, 10, [1, 2, 3]),
            Interval(1, 10, [1, 2, 3]),
            Interval(0, 10, [1, 2, 3])
        ]

        interval_collection = IntervalCollection(intervals)
        intervals_filtered = list(UniqueIntervals(interval_collection))

        self.assertEqual(len(intervals_filtered), len(intervals)-1)
        self.assertEqual(intervals_filtered[0], intervals[0])
        self.assertEqual(intervals_filtered[1], intervals[1])


if __name__ == "__main__":
    unittest.main()
