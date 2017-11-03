import unittest
import offsetbasedgraph as obg
from test_snarls import snarl_graph1, snarl_graph2
from graph_peak_caller.snarlmaps import LinearSnarlMap

graph = obg.GraphWithReversals(
    {3: obg.Block(20), 5: obg.Block(10),
     12: obg.Block(20), 13: obg.Block(21),
     }, {})


class TestSnarlMap(unittest.TestCase):
    def setUp(self):
        self.snarl_map = LinearSnarlMap(snarl_graph2, graph)
        self.graph_positions = [obg.Position(5, 4),
                                obg.Position(3, 4),
                                obg.Position(12, 4),
                                obg.Position(13, 4)]

        self.linear_positions = [4, 31/20*4, 10+21/20*4, 14]
        self.linear_positions = [int(p) for p in self.linear_positions]
        self.graph_interval = obg.DirectedInterval(self.graph_positions[0],
                                                   self.graph_positions[2])

    def test_graph_position_to_linear(self):
        for graph_pos, lin_pos in zip(self.graph_positions,
                                      self.linear_positions):
            mapped_pos = self.snarl_map.graph_position_to_linear(graph_pos)
            self.assertEqual(mapped_pos, lin_pos)

    def test_map_graph_interval(self):
        mapped_interval = self.snarl_map.map_graph_interval(
            self.graph_interval)
        self.assertEqual(mapped_interval, (self.linear_positions[0],
                                           self.linear_positions[2]), [5, 12])


if __name__ == "__main__":
    unittest.main()
