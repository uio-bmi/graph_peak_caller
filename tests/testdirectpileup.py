import numpy as np
from graph_peak_caller.samplepileup import PileupCreator, ReversePileupCreator
from graph_peak_caller.directsamplepileup import main, DirectPileup, Starts
from graph_peak_caller.samplepileup import sorted_wierd_graph as swierd
from graph_peak_caller.densepileup import DensePileup
import offsetbasedgraph as obg
# import cProfile
import unittest


def sorted_wierd_graph(a, b):
    nodes = {i+1: obg.Block(100) for i in range(4)}
    edges = {1: [2, 3], 2: [4], 3: [4]}
    return obg.GraphWithReversals(nodes, edges)


class TestDirectPileup(unittest.TestCase):
    def setUp(self):
        self._graph = sorted_wierd_graph(2, 100)
        self._pileup = DensePileup(self._graph)
        self._intervals = [obg.DirectedInterval(0, 10, [1]),
                           obg.DirectedInterval(90, 10, [1, 2]),
                           obg.DirectedInterval(90, 10, [-3, -1])]

    def test_direct(self):
        direct_pileup = DirectPileup(
            self._graph, self._intervals, self._pileup)
        direct_pileup.run()
        true_array = np.zeros(400, dtype="int")
        true_array[0:10] = 1
        true_array[90:100] = 2
        true_array[100:110] = 1
        true_array[200:210] = 1
        self.assertTrue(np.all(true_array == self._pileup.data._values))

    def test_end_postions(self):
        direct_pileup = DirectPileup(
            self._graph, self._intervals, self._pileup)
        direct_pileup.run()
        true_pos_ends = {1: [10], 2: [10], 3: [], 4: []}
        true_neg_ends = {-1: [10], -2: [], -3: [], -4: []}
        self.assertEqual(direct_pileup._pos_ends, true_pos_ends)
        self.assertEqual(direct_pileup._neg_ends, true_neg_ends)


class TestSamplePileup(unittest.TestCase):
    true_pos_ends = Starts({1: [10], 2: [10], 3: [], 4: []})
    true_neg_ends = Starts({-1: [10], -2: [], -3: [], -4: []})

    def setUp(self):
        self._graph = sorted_wierd_graph(2, 100)
        self._pileup = np.zeros(401, dtype="int")

    def test_pos_dir(self):
        creator = PileupCreator(
            self._graph, self.true_pos_ends, self._pileup)
        creator._fragment_length = 50
        creator.run_linear()
        true_array = np.zeros_like(creator._pileup)
        true_array[10:60] = 1
        true_array[110:160] = 1
        self.assertTrue(np.all(creator.get_pileup() == true_array))

    def test_neg_dir(self):
        creator = ReversePileupCreator(
            self._graph, self.true_neg_ends, self._pileup)
        creator._fragment_length = 50
        creator.run_linear()
        pileup = creator.get_pileup()
        true_array = np.zeros_like(pileup)
        true_array[40:90] = 1
        self.assertTrue(np.all(pileup == true_array))

    def test_pos_dir_long(self):
        creator = PileupCreator(
            self._graph, self.true_pos_ends, self._pileup)
        creator._fragment_length = 100
        creator.run_linear()
        true_array = np.zeros_like(creator._pileup)
        true_array[10:210] = 1
        true_array[300:310] = 1
        true_array[400:410] = 1
        print(creator.get_pileup() -true_array)
        self.assertTrue(np.all(creator.get_pileup()[:-1] == true_array[:-1]))

    def test_pos_dir_three_nodes(self):
        true_pos_ends = Starts({1: [90], 2: [70], 3: [], 4: []})
        creator = PileupCreator(
            self._graph, true_pos_ends, self._pileup)
        creator._fragment_length = 120
        creator.run_linear()
        true_array = np.zeros_like(creator._pileup)
        true_array[90:170] = 1
        true_array[170:200] = 2
        true_array[200:300] = 1
        true_array[300:310] = 2
        true_array[310:390] = 1
        print(creator.get_pileup()-true_array)
        self.assertTrue(np.all(creator.get_pileup() == true_array))

    def test_neg_dir_three_nodes(self):
        true_neg_ends = Starts({-4: [90], -3: [70], -2: [], -1: []})
        creator = ReversePileupCreator(
            self._graph, true_neg_ends, self._pileup)
        creator._fragment_length = 120
        creator.run_linear()
        pileup = creator.get_pileup()
        true_array = np.zeros_like(pileup)
        true_array[300:310] = 1
        true_array[230:300] = 1
        true_array[200:230] = 2
        true_array[100:200] = 1
        true_array[90:100] = 2
        true_array[10:90] = 1
        self.assertTrue(np.all(pileup == true_array))


class TestMain(unittest.TestCase):
    def setUp(self):
        self._graph = sorted_wierd_graph(2, 100)
        self._pileup = DensePileup(self._graph)
        self._intervals = [obg.DirectedInterval(0, 20, [1]),
                           obg.DirectedInterval(90, 10, [1, 2]),
                           obg.DirectedInterval(90, 10, [-3, -1])]

    def test_simple(self):
        pileup = main(self._intervals, self._graph, 30)
        pileup = pileup.data._values
        true_array = np.zeros_like(pileup)
        true_array[0:50] = 1
        true_array[60:90] = 1
        true_array[90:100] = 2
        true_array[100:140] = 1
        true_array[200:210] = 1
        self.assertTrue(np.all(pileup == true_array))

    def test_long_fragments(self):
        intervals = [obg.DirectedInterval(80, 90, [1]),
                     obg.DirectedInterval(90, 10, [1, 2]),
                     obg.DirectedInterval(90, 10, [-4, -2]),
                     obg.DirectedInterval(90, 10, [-3, -1])]
        pileup = main(intervals, self._graph, 30)
        pileup = pileup.data._values
        true_array = np.zeros_like(pileup)
        true_array[60:80] = 1
        true_array[80:90] = 2
        true_array[90:100] = 3

        true_array[100:120] = 2
        true_array[120:140] = 1
        true_array[160:200] = 1

        true_array[200:210] = 2
        true_array[210:220] = 1

        true_array[300:310] = 1
        self.assertTrue(np.all(pileup == true_array))


if False or __name__ == "__main__":
    node_size = 32
    n_nodes = 100000
    graph = swierd(n_nodes, node_size)
    intervals = []
    for node in graph.get_sorted_node_ids():
        if node % 20 == 0:
            next_nodes = graph.adj_list[node]
            if not next_nodes:
                continue
            next_node = next_nodes[0]
            for i in range(10):
                intervals.append(
                    obg.DirectedInterval(i, i+4, [node, next_node]))

    for node in graph.get_sorted_node_ids(reverse=True):
        node = -node
        if node % 20 == 0:
            next_nodes = graph.reverse_adj_list[node]
            if not next_nodes:
                continue
            next_node = next_nodes[0]
            for i in range(10):
                intervals.append(
                    obg.DirectedInterval(i, i+4, [node, next_node]))

    main(intervals, graph, 110)
