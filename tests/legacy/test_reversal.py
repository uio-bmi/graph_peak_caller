import unittest
import logging
import pytest
if pytest.__version__ < "3.0.0":
    pytest.skip()

# from graph_peak_caller.extender import Extender, Areas
from graph_peak_caller.areas import BinaryContinousAreas
from offsetbasedgraph import Block, GraphWithReversals,\
    DirectedInterval, Position
from cyclic_graph import get_small_cyclic_graph, get_large_cyclic_graph
logging.getLogger("extender").setLevel("DEBUG")


logging.basicConfig(level=logging.DEBUG)

Graph = GraphWithReversals
nodes = {i: Block(20) for i in range(1, 4)}
tmp_edges = {1: [2, -2],
             2: [3],
             -2: [3]}
graph = Graph(nodes, tmp_edges)

pos_interval = DirectedInterval(8, 18, [2])
neg_interval = DirectedInterval(2, 12, [-2])
neg_spanning_interval = DirectedInterval(12, 2, [-3, -2])
pos_spanning_interval = DirectedInterval(18, 8, [2, 3])

@pytest.mark.skip("Legacy")
class TestExtender(unittest.TestCase):

    def setUp(self):
        self.extender = Extender(graph, 20)

    def test_extend(self):
        areas = self.extender.extend_interval(pos_interval)
        true_areas = BinaryContinousAreas(graph)
        true_areas.add_start(-2, 12)
        true_areas.add_start(3, 8)
        true_areas.add_start(-1, 8)
        self.assertEqual(areas, true_areas)

    def test_extend_reverse(self):
        areas = self.extender.extend_interval(neg_interval)
        true_areas = BinaryContinousAreas(graph)
        true_areas.add_start(2, 18)
        true_areas.add_start(3, 2)
        true_areas.add_start(-1, 2)
        self.assertEqual(areas, true_areas)

    def test_extend_cyclic(self):
        graph = get_small_cyclic_graph()
        interval = DirectedInterval(70, 90, [1], graph=graph)
        extender = Extender(graph, 40)
        areas = extender.extend_interval(interval)
        true_areas = BinaryContinousAreas(graph)
        true_areas.add_start(1, 10)
        true_areas.add_start(-1, 30)

        self.assertEqual(areas, true_areas)

    def test_extend_large_cyclic(self):
        logging.debug("CYCLIC")
        graph = get_large_cyclic_graph()
        interval = DirectedInterval(90, 10, [1, 2], graph=graph)
        extender = Extender(graph, 40)
        areas = extender.extend_interval(interval)
        true_areas = BinaryContinousAreas(graph)
        true_areas.add_start(1, 10)
        true_areas.add_start(-1, 10)
        true_areas.add_start(2, 20)
        self.assertEqual(areas, true_areas)

    def _test_areas_from_point_pos(self):
        traverser = self.extender.pos_traverser
        point = Position(2, 8)
        areas = self.extender.get_areas_from_point(point, 20, traverser)
        true_areas = Areas(graph, {2: [8, 20],
                                   3: [0, 8]})
        self.assertEqual(areas, true_areas)

    def _test_areas_from_point_neg(self):
        traverser = self.extender.neg_traverser
        point = Position(-2, 12)
        areas = self.extender.get_areas_from_point(point, 20, traverser)
        true_areas = Areas(graph, {-2: [12, 20],
                                   -1: [0, 12]})
        self.assertEqual(areas, true_areas)


@pytest.mark.skip("Legacy")
class TestAreas(unittest.TestCase):
    def test_from_pos_interval(self):
        areas = Areas.from_interval(pos_interval, graph)
        true_areas = Areas(graph, {2: [8, 18]})
        self.assertEqual(areas, true_areas)

    def test_from_neg_interval(self):
        areas = Areas.from_interval(neg_interval, graph)
        areas.reverse_reversals()
        true_areas = Areas(graph, {2: [8, 18]})
        self.assertEqual(areas, true_areas)

    def test_from_spanning_neg_interval(self):
        areas = Areas.from_interval(neg_spanning_interval, graph)
        areas.reverse_reversals()
        true_areas = Areas(graph, {2: [18, 20],
                                   3: [0, 8]})
        self.assertEqual(areas, true_areas)

    def test_from_spanning_pos_interval(self):
        areas = Areas.from_interval(pos_spanning_interval, graph)
        areas.reverse_reversals()
        true_areas = Areas(graph, {2: [18, 20],
                                   3: [0, 8]})
        self.assertEqual(areas, true_areas)

    def test_update_areas(self):
        areas1 = Areas(graph, {2: [0, 20]})
        areas2 = Areas(graph, {2: [0, 10]})
        areas1.update(areas2)
        self.assertEqual(areas1, Areas(graph, {2: [0, 20]}))

    def test_update_areas_mid(self):
        areas1 = Areas(graph, {2: [10, 20]})
        areas2 = Areas(graph, {2: [0, 10]})
        areas1.update(areas2)
        self.assertEqual(areas1, Areas(graph, {2: [0, 20]}))

    def test_update_areas_mid2(self):
        areas1 = Areas(graph, {2: [5, 12]})
        areas2 = Areas(graph, {2: [12, 17]})
        areas1.update(areas2)
        self.assertEqual(areas1, Areas(graph, {2: [5, 17]}))

    def test_reverse_reversals(self):
        areas = Areas(graph, {2: [8, 18],
                              -2: [2, 12]})
        areas.reverse_reversals()
        self.assertEqual(areas, Areas(graph, {2: [8, 18]}))


if __name__ == "__main__":
    unittest.main()
