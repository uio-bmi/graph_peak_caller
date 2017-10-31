import unittest
from offsetbasedgraph import Graph, Block
from graph_peak_caller.snarls import SimpleSnarl, SnarlGraphBuilder, SnarlGraph

class SnarlsTests(unittest.TestCase):
    def setUp(self):
        self.simple_graph = Graph(
            {i: Block(3) for i in range(1, 9)},
            {
                1: [2, 3],
                 2: [4],
                 3: [4],
                 4: [5],
                 5: [6, 7],
                 6: [8],
                 7: [8]
             })

        self.simple_snarls = \
            {
                10: SimpleSnarl(1, 4, 10),
                11: SimpleSnarl(5, 8, 11)
            }

class TestSnarlGraphBuilder(SnarlsTests):

    def test_build_non_nested(self):

        builder = SnarlGraphBuilder(self.simple_graph, self.simple_snarls)
        graph = builder.build_snarl_graphs()

        print(graph)





if __name__ == "__main__":
    unittest.main()
