import unittest

from graph_peak_caller.areas import BinaryContinousAreas, BCACollection
import offsetbasedgraph as obg

nodes = {i: obg.Block(10) for i in range(1, 10)}
graph = obg.GraphWithReversals(nodes, {})


class TestBinaryContinousAreas(unittest.TestCase):

    def setUp(self):
        self.areas = BinaryContinousAreas(graph)
        self.areas.full_areas = {1: 1, 2: 1}
        self.areas.starts = {3: 2, 4: 5}
        self.areas.internal_intervals = {5: [2, 7]}

    def test_file_in_out(self):
        line = self.areas.to_file_line()
        print(line)
        new_areas = BinaryContinousAreas.from_file_line(line, graph)
        self.assertEqual(new_areas, self.areas)
        c = BCACollection([self.areas, self.areas])
        c.to_file("tmp.subgraphs")
        BCACollection.from_file("tmp.subgraphs", graph)

if __name__ == "__main__":
    unittest.main()
