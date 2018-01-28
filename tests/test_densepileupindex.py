import unittest
from graph_peak_caller.densepileupindex import DensePileupExtender, GraphIndex, GraphExtender
from offsetbasedgraph import GraphWithReversals as Graph, Block, DirectedInterval as Interval
from graph_peak_caller.densepileup import DensePileup
from graph_peak_caller.sampleandcontrolcreator import create_sample_using_indexes


class TestGraphExtenderLinearGraph(unittest.TestCase):

    def setUp(self):
        self.graph = Graph({i: Block(10) for i in range(1, 4)},
                           {i: [i+1] for i in range(1, 3)})

        self.index = GraphIndex(
            {
                1: [(2, 10), (3, 20)],
                2: [(3, 10)],
                3: [],
                -1: [],
                -2: [(-1, 10)],
                -3: [(-2, 10), (-1, 20)]
            }
        )
        self.extender = GraphExtender(self.index)

    def test_simple(self):
        extensions = self.extender.extend_from_position(1, 4, extension_length=1)
        self.assertEqual(list(extensions),  [(1, -4)])

        extensions = self.extender.extend_from_position(1, 3, extension_length=1)
        self.assertEqual(list(extensions),  [(1, -3)])

        extensions = self.extender.extend_from_position(1, 3, extension_length=10)
        self.assertEqual(list(extensions),  [(1, -3), (2, 7)])

        extensions = self.extender.extend_from_position(1, 3, extension_length=20)
        self.assertEqual(list(extensions),  [(1, -3), (2, 7), (3, 17)])

        extensions = self.extender.extend_from_position(1, 3, extension_length=100)
        self.assertEqual(list(extensions),  [(1, -3), (2, 7), (3, 17)])

    def test_reverse_extension(self):
        extensions = self.extender.extend_from_position(-1, 4, extension_length=1)
        self.assertEqual(list(extensions),  [(-1, -4)])

        extensions = self.extender.extend_from_position(-3, 2, extension_length=20)
        self.assertEqual(list(extensions),  [(-3, -2), (-2, 8), (-1, 18)])


class TestGraphExtenderSplitGraph(unittest.TestCase):

    def setUp(self):
        self.graph = Graph({i: Block(10) for i in range(1, 5)},
                           {
                               1: [2, 3],
                               2: [4],
                               3: [4]
                           })

        self.index = GraphIndex(
            {
                1: [(2, 10), (3, 10), (4, 20)],
                2: [(4, 10)],
                3: [(4, 10)],
                4: [],
                -1: [],
                -2: [(-1, 10)],
                -3: [(-1, 10)],
                -4: [(-2, 10), (-3, 10), (-1, 20)]
            }
        )
        self.extender = GraphExtender(self.index)

    def test_forward(self):
        extensions = self.extender.extend_from_position(1, 4, extension_length=10)
        self.assertEqual(list(extensions),  [(1, -4), (2, 6), (3, 6)])

        extensions = self.extender.extend_from_position(1, 5, extension_length=20)
        self.assertEqual(list(extensions),  [(1, -5), (2, 5), (3, 5), (4, 15)])


class TestPileupExtender(unittest.TestCase):

    def setUp(self):
        self.graph = Graph({i: Block(10) for i in range(1, 4)},
                           {i: [i+1] for i in range(1, 3)})

        self.index = GraphIndex(
            {
                1: [(2, 10), (3, 20)],
                2: [(3, 10)],
                3: [],
                -1: [],
                -2: [(-1, 10)],
                -3: [(-2, 10), (-1, 20)]
            }
        )
        pileup = DensePileup(self.graph)
        self.extender = DensePileupExtender(GraphExtender(self.index), pileup)

    def test_simple(self):
        extensions = self.extender.extend_from_position(1, 4, extension_length=1)
        self.assertEqual(list(extensions),  [(4, 5)])

        extensions = self.extender.extend_from_position(1, 4, extension_length=2)
        self.assertEqual(list(extensions),  [(4, 6)])

        extensions = self.extender.extend_from_position(1, 4, extension_length=6)
        self.assertEqual(list(extensions),  [(4, 10)])

        extensions = self.extender.extend_from_position(1, 4, extension_length=7)
        self.assertEqual(list(extensions),  [(4, 10), (10, 11)])

        extensions = self.extender.extend_from_position(1, 4, extension_length=8)
        self.assertEqual(list(extensions),  [(4, 10), (10, 12)])

    def test_reverse(self):
        extensions = self.extender.extend_from_position(-1, 3, extension_length=2)
        self.assertEqual(list(extensions),  [(5, 7)])

        extensions = self.extender.extend_from_position(-1, 3, extension_length=7)
        self.assertEqual(list(extensions),  [(0, 7)])

        extensions = self.extender.extend_from_position(-2, 3, extension_length=17)
        self.assertEqual(list(extensions),  [(10, 17), (0, 10)])


class TestGraphIndex(unittest.TestCase):
    def test_from_linear_graph(self):
        graph = Graph({i: Block(10) for i in range(1, 4)},
                           {i: [i+1] for i in range(1, 3)})

        index = GraphIndex.create_from_graph(graph, length=25)
        self.assertEqual(list(index.get_node(1)), [(1, 0), (2, 10), (3, 20)])
        self.assertEqual(list(index.get_node(2)), [(2, 0), (3, 10)])
        self.assertEqual(list(index.get_node(3)), [(3, 0)])
        self.assertEqual(list(index.get_node(-1)), [(-1, 0)])
        self.assertEqual(list(index.get_node(-2)), [(-2, 0), (-1, 10)])
        self.assertEqual(list(index.get_node(-3)), [(-3, 0), (-1, 20), (-2, 10)])

    def test_from_split_graph(self):
        graph = Graph(
            {1: Block(10),
             2: Block(3),
             3: Block(8),
             4: Block(10),
             },
            {
            1: [2, 3],
            2: [4],
            3: [4]
            }
        )
        index = GraphIndex.create_from_graph(graph, length=15)
        self.assertEqual(list(index.get_node(1)), [(1, 0), (2, 10), (3, 10), (4, 13)])

        index = GraphIndex.create_from_graph(graph, length=11)
        self.assertEqual(list(index.get_node(1)), [(1, 0), (2, 10), (3, 10)])

    def test_to_from_file(self):
        index = GraphIndex(
            {
                1: [(2, 123), (3, 4)],
                2: [(3, 41)]
            }
        )
        index.to_file("test")
        read_index = GraphIndex.from_file("test")
        self.assertEqual(index, read_index)


class TestCreateSamplePileupUsingGraphIndex(unittest.TestCase):

    def setUp(self):

        self.linear_graph = Graph({i: Block(10) for i in range(1, 4)},
                           {i: [i+1] for i in range(1, 3)})
        self.linear_graphindex = GraphIndex.create_from_graph(self.linear_graph, 100)

        self.split_graph = Graph(
            {1: Block(10),
             2: Block(3),
             3: Block(8),
             4: Block(10),
             },
            {
            1: [2, 3],
            2: [4],
            3: [4]
            }
        )
        self.split_graphindex = GraphIndex.create_from_graph(self.split_graph, length=100)


    def test_linear_graph(self):
        graph = self.linear_graph
        graphindex = self.linear_graphindex
        sample_reads = [
            Interval(5, 7, [1], graph)
        ]
        pileup = create_sample_using_indexes(sample_reads, graphindex, graph, fragment_length=5)
        self.assertEqual(pileup, DensePileup.from_intervals(graph, [Interval(5, 10, [1])]))

        sample_reads = [
            Interval(5, 7, [1], graph)
        ]
        pileup = create_sample_using_indexes(sample_reads, graphindex, graph, fragment_length=5)
        self.assertEqual(pileup, DensePileup.from_intervals(graph, [Interval(5, 10, [1])]))

        sample_reads = [
            Interval(3, 5, [-1], graph)
        ]
        pileup = create_sample_using_indexes(sample_reads, graphindex, graph, fragment_length=5)
        self.assertEqual(pileup, DensePileup.from_intervals(graph, [Interval(2, 7, [1])]))

    def test_linear_graph_reverse(self):
        graph = self.linear_graph
        graphindex = self.linear_graphindex
        sample_reads = [
            Interval(8, 2, [-3, -2], graph)
        ]
        pileup = create_sample_using_indexes(sample_reads, graphindex, graph, fragment_length=15)
        self.assertEqual(pileup, DensePileup.from_intervals(graph, [Interval(7, 2, [1, 2, 3])]))

    def test_split_graph(self):
        graph = self.split_graph
        graphindex = self.split_graphindex

        sample_reads = [
            Interval(8, 9, [1], graph)
        ]
        pileup = create_sample_using_indexes(sample_reads, graphindex, graph, fragment_length=15)
        self.assertEqual(pileup, DensePileup.from_intervals(graph, [Interval(8, 10, [1, 2, 4]), Interval(0, 8, [3])]))


if __name__ == "__main__":
    unittest.main()
