import unittest
from graph_peak_caller.haplotyping import HaploTyper
from offsetbasedgraph import DirectedInterval as Interval, GraphWithReversals as Graph, Block, IntervalCollection

class TestHaplotyper(unittest.TestCase):

    def setUp(self):
        self.complex_graph = Graph(
            {i: Block(3) for i in range(1, 13)},
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
                10: [12]
             })
        self.complex_graph.convert_to_numpy_backend()

    def test_simple(self):
        graph = Graph(
            {i: Block(3) for i in range(1, 5)},
            {
                1: [2, 3],
                2: [4],
                3: [4]
            }
        )
        graph.convert_to_numpy_backend()

        intervals = IntervalCollection([
            Interval(0, 3, [1, 3])
        ])

        haplotyper = HaploTyper(graph, intervals)
        haplotyper.build()
        max_interval = haplotyper.get_maximum_interval_through_graph()

        self.assertEqual(
            max_interval,
            Interval(0, 3, [1, 3, 4])
        )

    def test_complex_graph(self):
        intervals = IntervalCollection([
            Interval(0, 3, [1, 3, 4, 6, 10]),
            Interval(1, 2, [2]),
            Interval(2, 3, [2]),
            Interval(0, 3, [7, 9])
        ])
        haplotyper = HaploTyper(self.complex_graph, intervals)
        haplotyper.build()
        max_interval = haplotyper.get_maximum_interval_through_graph()

        self.assertEqual(
            max_interval,
            Interval(0, 3, [1, 2, 7, 9, 10, 12])
        )

if __name__ == "__main__":
    unittest.main()