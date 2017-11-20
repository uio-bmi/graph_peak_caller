import unittest

from graph_peak_caller.areas import BinaryContinousAreas, BCACollection
import offsetbasedgraph as obg

nodes = {i: obg.Block(10) for i in range(1, 10)}
graph = obg.GraphWithReversals(nodes, {i: [i+1] for i in range(1, 9)})


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

    def test_filled_interval(self):
        interval = obg.DirectedInterval(4, 4, [2, 3, 4])
        areas = BinaryContinousAreas(graph)
        areas.filled_interval(interval)
        areas.sanitize()
        true_areas = BinaryContinousAreas(graph)
        true_areas.full_areas = {3: 1}
        true_areas.starts = {-2: 6, 4: 4}
        self.assertEqual(areas, true_areas)

    def test_internal_filled_interval(self):
        interval = obg.DirectedInterval(2, 8, [3])
        areas = BinaryContinousAreas(graph)
        areas.filled_interval(interval)
        areas.sanitize()
        true_areas = BinaryContinousAreas(graph)
        true_areas.internal_intervals = {3: [2, 8]}
        self.assertEqual(areas, true_areas)


if __name__ == "__main__":
    unittest.main()
