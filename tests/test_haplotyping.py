import unittest
from graph_peak_caller.haplotyping import HaploTyper
from offsetbasedgraph import DirectedInterval, GraphWithReversals as Graph, Block, IntervalCollection

class TestHaplotyper(unittest.TestCase):

    def test_simple(self):
        graph = Graph(
            {i: Block(3) for i in range(1, 5)},
            {
                1: [2, 3],
                2: [4],
                3: [4]
            }
        )

        intervals = IntervalCollection([
            DirectedInterval(0, 3, [1, 3])
        ])

        haplotyper = HaploTyper(graph, intervals)
        haplotyper.build()
        max_interval = haplotyper.get_maximum_interval_through_graph()

        self.assertEqual(
            max_interval,
            DirectedInterval(0, 3, [1, 3, 4])
        )

if __name__ == "__main__":
    unittest.main()