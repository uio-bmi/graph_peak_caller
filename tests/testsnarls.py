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
                20: SimpleSnarl(1, 4, 20),
                21: SimpleSnarl(5, 8, 21),
                22: SimpleSnarl(4, 5, 22)
            }

class TestSnarlGraphBuilder(SnarlsTests):

    def test_build_non_nested(self):

        builder = SnarlGraphBuilder(self.simple_graph, self.simple_snarls)
        graph = builder.build_snarl_graphs()

        correct_snarl_graph = SnarlGraph(
            {1: Block(3),
             20: SnarlGraph(
                 {
                     2: Block(3),
                     3: Block(3)
                 },
                 {}
             ),
             4:  Block(3),
             5:  Block(3),
             21: SnarlGraph(
                 {
                     6: Block(3),
                     7: Block(3)
                 },
                 {}
             ),
             8: Block(3)
             },
            {
                1: [20],
                20: [4],
                4: [5],
                5: [21],
                21: [8],
                9: [1],  # Dummy nodes added by snarlgraphbuilder
                8: [10]   # Dummy nodes added by snarlgraphbuilder
            }
        )
        print("Final graph")

        print(graph)

        self.assertEqual(correct_snarl_graph, graph)


    def test_hierarchical(self):
        graph = Graph(
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

        subsnarl1 = SimpleSnarl(3, 6, 21, parent=20)
        subsnarl2 = SimpleSnarl(2, 9, 22, parent=20)
        parent_snarl = SimpleSnarl(1, 10, 20, children=[subsnarl1, subsnarl2])

        snarls = {
            20: parent_snarl,
            21: subsnarl1,
            22: subsnarl2
        }

        builder = SnarlGraphBuilder(graph, snarls)
        snarlgraph = builder.build_snarl_graphs()
        print("Snarlgraph")
        print(snarlgraph)

        correct_snarl_graph = SnarlGraph(
            {
                11: Block(3),
                12: Block(3),
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
                            {}
                        ),
                        22: SnarlGraph(
                            {
                                7: Block(3),
                                8: Block(3)
                            },
                            {}
                        ),
                        2: Block(3),
                        6: Block(3),
                        9: Block(3),
                    },
                    {
                        3: [21],
                        2: [22],
                        21: [6],
                        22: [9]
                    }
                )
            },
            {
                11: [1],
                1: [20],
                20: [10],
                10: [12],
                13: [11],  # Dummy
                12: [14]   # Dummy
            }
        )

        print("Snarlgraph")
        print(snarlgraph)

        self.assertEqual(correct_snarl_graph, snarlgraph)

    def _test_multiple_snarls_same_nodes(self):

        graph = Graph(
            {i: Block(3) for i in range(1, 7)},
            {
                1: [2, 4],
                2: [3],
                3: [6],
                4: [5],
                5: [6]
            })

        snarls = {
            10: SimpleSnarl(1, 3, 10),
            11: SimpleSnarl(1, 5, 11),
        }

        builder = SnarlGraphBuilder(graph, snarls)
        snarlgraph = builder.build_snarl_graphs()

        correct_snarl_graph = SnarlGraph(
            {
                1: Block(3),
                10: SnarlGraph(
                    {
                        2: Block(3)
                    }, {}
                ),
                11: SnarlGraph(
                    {
                        4: Block(3)
                    }, {}
                ),
                4: Block(3),
                5: Block(3),
            },
            {
                1: [10, 11],
                10: [3],
                11: [5],
                3: [6],
                5: [6]
            })

        self.assertEqual(correct_snarl_graph, snarlgraph)





if __name__ == "__main__":
    unittest.main()
