import unittest
from graph_peak_caller.densepileupindex import DensePileupExtender, GraphIndex, GraphExtender
from offsetbasedgraph import GraphWithReversals as Graph, Block, DirectedInterval as Interval
from graph_peak_caller.densepileup import DensePileup


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

if __name__ == "__main__":
    unittest.main()
