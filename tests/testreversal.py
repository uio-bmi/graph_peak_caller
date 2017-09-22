import unittest
from collections import defaultdict
from graph_peak_caller.extender import Extender, Areas
from offsetbasedgraph import Block, Graph, Interval
nodes = {i: Block(20) for i in range(4)}
tmp_edges = {1: [2, -2],
             2: [3],
             -2: [3]}

edges = defaultdict(list)
edges.update(tmp_edges)

tmp_rev_edges = {-2: [-1],
                 2: [-1],
                 3: [-2, 2]}
rev_edges = defaultdict(list)
rev_edges.update(tmp_rev_edges)
graph = Graph(edges, tmp_edges)

pos_interval = Interval(5, 15, [2])


class TestExtender(unittest.TestCase):
    def test_extend(self):
        extender = Extender(graph, 10)
        areas = extender.extend_interval(interval)
        true_areas = Areas(graph, {2: [5, 20],
                                   3: [0, 5]})
        self.assertEqual(areas, true_areas)

