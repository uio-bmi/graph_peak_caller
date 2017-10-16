import unittest
from graph_peak_caller.subgraphcollection import SubgraphCollection, ConnectedAreas
from graph_peak_caller.extender import Areas
from offsetbasedgraph import Graph, Block, Interval, Position
import numpy as np
from graph_peak_caller.sparsepileup import SparsePileup

class Tester(unittest.TestCase):
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

        areas = {
            2: np.array([0, 3])
        }
        self.middle_areas = ConnectedAreas(self.simple_graph, areas)
        self.middle_closed_area = ConnectedAreas(self.simple_graph, {2: np.array([1, 2])})
        self.middle_left_area = ConnectedAreas(self.simple_graph, {2: np.array([0, 2])})

class TestSubGraphCollection(Tester):

    def test_add_single_area(self):
        collection = SubgraphCollection(self.simple_graph)
        collection.add_area(1, 2, 3)

        self.assertEqual(len(collection.subgraphs), 1)
        self.assertTrue(1 in collection.subgraphs[0].areas)
        self.assertTrue(np.all(np.array([2, 3]) == collection.subgraphs[0].areas[1]))

    def test_add_two_disjoint_areas(self):
        collection = SubgraphCollection(self.simple_graph)
        collection.add_area(1, 1, 2)
        collection.add_area(2, 2, 3)

        self.assertEqual(len(collection.subgraphs), 2)
        self.assertTrue(1 in collection.subgraphs[0].areas)
        self.assertTrue(2 in collection.subgraphs[1].areas)
        self.assertTrue(np.all(np.array([1, 2]) == collection.subgraphs[0].areas[1]))
        self.assertTrue(np.all(np.array([2, 3]) == collection.subgraphs[1].areas[2]))

    def test_add_two_connected_areas(self):
        print(self.simple_graph.adj_list)
        print(self.simple_graph.reverse_adj_list)
        collection = SubgraphCollection(self.simple_graph)
        collection.add_area(1, 2, 3)
        collection.add_area(2, 0, 3)

        self.assertTrue(len(collection.subgraphs), 1)
        areas = collection.subgraphs[0].areas
        self.assertTrue(np.all(areas[1] == [2, 3]))
        self.assertTrue(np.all(areas[2] == [0, 3]))

    def test_two_subgraphs_where_one_is_connected(self):
        collection = SubgraphCollection(self.simple_graph)
        collection.add_area(1, 2, 3)
        collection.add_area(2, 0, 2)
        collection.add_area(3, 1, 2)

        self.assertTrue(len(collection.subgraphs), 2)
        areas = collection.subgraphs[0].areas
        self.assertTrue(np.all(areas[1] == [2, 3]))
        self.assertTrue(np.all(areas[2] == [0, 2]))

        areas = collection.subgraphs[1].areas
        self.assertTrue(np.all(areas[3] == [1, 2]))

    def test_subgraphs_with_loop(self):
        collection = SubgraphCollection(self.simple_graph)
        collection.add_area(1, 2, 3)
        collection.add_area(2, 0, 2)
        collection.add_area(3, 1, 2)
        collection.add_area(1, 0, 1)

        self.assertTrue(len(collection.subgraphs), 3)
        areas = collection.subgraphs[0].areas
        print(areas[1])
        self.assertTrue(np.all(areas[1] == [2, 3]))
        self.assertTrue(np.all(areas[2] == [0, 2]))

        areas = collection.subgraphs[1].areas
        self.assertTrue(np.all(areas[3] == [1, 2]))

        areas = collection.subgraphs[2].areas
        self.assertTrue(np.all(areas[1] == [0, 1]))

    def test_subgraphs_multiple(self):
        collection = SubgraphCollection(self.simple_graph)
        collection.add_area(1, 0, 3)
        collection.add_area(3, 0, 3)
        collection.add_area(2, 0, 3)

        self.assertTrue(len(collection.subgraphs), 1)
        areas = collection.subgraphs[0].areas
        self.assertTrue(np.all(areas[1] == [0, 3]))
        self.assertTrue(np.all(areas[2] == [0, 3]))
        self.assertTrue(np.all(areas[3] == [0, 3]))

    def test_create_subgraph_from_pileup(self):

        intervals = [Interval(0, 3, [1, 2, 3], self.simple_graph)]
        pileup = SparsePileup.from_intervals(self.simple_graph, intervals)
        collection = SubgraphCollection.from_pileup(self.simple_graph, pileup)
        self.assertTrue(len(collection.subgraphs), 1)
        areas = collection.subgraphs[0].areas
        self.assertTrue(np.all(areas[1] == [0, 3]))
        self.assertTrue(np.all(areas[2] == [0, 3]))
        self.assertTrue(np.all(areas[3] == [0, 3]))



class TestConnectedAreas(Tester):

    def test_connected_area_touches(self):

        self.assertTrue(self.middle_areas.touches_area(1, 0, 3))
        self.assertTrue(self.middle_areas.touches_area(3, 0, 3))
        self.assertTrue(self.middle_areas.touches_area(3, 0, 2))

        self.assertFalse(self.middle_closed_area.touches_area(3, 0, 3))
        self.assertFalse(self.middle_closed_area.touches_area(1, 0, 3))

        self.assertFalse(self.middle_left_area.touches_area(3, 0, 3))
        self.assertTrue(self.middle_left_area.touches_area(1, 0, 3))

if __name__ == "__main__":
    unittest.main()