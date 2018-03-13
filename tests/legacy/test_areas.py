import unittest
import pytest
if pytest.__version__ < "3.0.0":
    pytest.skip()

# from graph_peak_caller.subgraphcollection import SubgraphCollection, ConnectedAreas
# from graph_peak_caller.extender import Areas
# from offsetbasedgraph import Graph, Block, Interval, Position
# import numpy as np
# from graph_peak_caller.legacy.sparsepileup import SparsePileup


@pytest.mark.skip("Legacy")
class TestAreas(unittest.TestCase):

    def setUp(self):
        self.simple_graph = Graph({
            1: Block(3),
            2: Block(3),
            3: Block(3)
        },
        {
            1: [2],
            2: [3]
        })

        self.reversed_graph = Graph({
            1: Block(3),
            2: Block(3),
            3: Block(3)
        },
        {
            -2: [-1],
            -3: [-2]
        })

    def test_get_start_and_end_positions(self):
        for graph in [self.simple_graph, self.reversed_graph]:
            areas = Areas(graph,
                          {1: np.array([2, 3]),
                           2: np.array([0, 3]),
                           3: np.array([0, 2])
                           }
                    )

            starts_ends = areas.get_start_and_end_positions()
            self.assertTrue(Position(1, 2) in starts_ends)
            self.assertTrue(Position(3, 1) in starts_ends)



if __name__ == "__main__":
    unittest.main()
