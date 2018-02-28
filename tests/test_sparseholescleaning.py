import numpy as np
import offsetbasedgraph as obg
from graph_peak_caller.sparseholecleaner import HolesCleaner
from graph_peak_caller.sparsediffs import SparseValues
import pytest


@pytest.fixture
def complicated_graph():
    nodes = {i: obg.Block(2) for i in range(1, 11)}
    edges = {1: [2, 3],
             2: [4],
             3: [4],
             4: [5, 6],
             5: [7],
             6: [7],
             7: [8, 9],
             9: [10]}
    return obg.GraphWithReversals(nodes, edges)


@pytest.fixture
def complicated_offset():
    nodes = {i: obg.Block(2) for i in range(101, 111)}
    edges = {101: [102, 103],
             102: [104],
             103: [104],
             104: [105, 106],
             105: [107],
             106: [107],
             107: [108, 109],
             109: [110]}
    return obg.GraphWithReversals(nodes, edges)


@pytest.fixture
def small_graph():
    nodes = {i: obg.Block(10) for i in range(101, 107)}
    edges = {101: [102],
             102: [103, 104],
             103: [105],
             104: [105],
             105: [106]}
    return obg.GraphWithReversals(nodes, edges)


def test_holes_cleaner():
    indices = np.array([80, 100,
                        180, 220,
                        240, 250,
                        300, 400,
                        500, 520,
                        610, 810])
    values = np.array([(i % 2) for i, _ in enumerate(indices)])
    pileup = SparseValues(indices, values)
    graph = obg.GraphWithReversals({i+1: obg.Block(100) for i in range(10)},
                                   {i: [i+1] for i in range(1, 10)})
    # graph.node_indexes = np.arange(0, 1001, 100)
    holes = HolesCleaner(graph, pileup, 10).run()
    print(holes)


def test_end_hole():
    pileup = SparseValues([0, 5, 10],
                          [False, True, False])
    graph = obg.GraphWithReversals(
        {1: obg.Block(12)}, {1: []})
    holes = HolesCleaner(graph, pileup, 4).run()
    assert holes == pileup


def test_long_hole():
    pileup = SparseValues([0, 5, 95],
                          [True, False, True])
    graph = obg.GraphWithReversals(
        {i: obg.Block(1) for i in range(1, 101)},
        {i: [i+1] for i in range(1, 100)})
    holes = HolesCleaner(graph, pileup, 56).run()
    assert holes == pileup


def test_complicated(complicated_graph):
    # graph = complicated_graph()
    #                      
    pileup = SparseValues(
        #         1  1  2--3  4  5---6   7   8
        np.array([0, 1, 2, 4, 6, 8, 10, 12, 14], dtype="int"),
        np.array([1, 0, 1, 0, 1, 0,  1, 0,  1], dtype="bool"))
    holes = HolesCleaner(complicated_graph, pileup, 3).run()
    true = SparseValues([0, 8, 10], [1, 0, 1])
    assert holes == true


def test_offset(complicated_offset):
    pileup = SparseValues(
        #         1  1  2--3  4  5---6   7   8
        np.array([0, 1, 2, 4, 6, 8, 10, 12, 14], dtype="int"),
        np.array([1, 0, 1, 0, 1, 0,  1, 0,  1], dtype="bool"))
    holes = HolesCleaner(complicated_offset, pileup, 3).run()
    true = SparseValues([0, 8, 10], [1, 0, 1])
    assert holes == true


def test_internals():
    graph = obg.GraphWithReversals({101: obg.Block(100)},
                                   {101: []})

    pileup = SparseValues(
        [0, 10, 19, 30, 41, 50],
        [1,  0,  1,  0,  1,  0])
    cleaned = HolesCleaner(graph, pileup, 10).run()
    true = SparseValues(
        [0, 30, 41, 50],
        [1, 0,  1,  0])

    assert cleaned == true


def test_touched_edge(small_graph):
    touched_nodes = set([101, 103, 106])
    pileup = SparseValues([0, 4, 22, 28, 56],
                          [1, 0,  1,  0,  1])
    cleaned = HolesCleaner(small_graph, pileup, 4, touched_nodes).run()
    assert cleaned == pileup


# @pytest.mark.skip
def test_non_touched_mid(small_graph):
    touched_nodes = set([101, 102, 103, 104, 105, 106])
    pileup = SparseValues([0, 4, 22, 28, 56],
                          [1, 0,  1,  0,  1])
    cleaned = HolesCleaner(small_graph, pileup, 20, touched_nodes).run()
    true = SparseValues([0, 30, 40], [1, 0, 1])
    assert cleaned == true


if __name__ == "__main__":
    # test_holes_cleaner()
    # test_end_hole()
    test_internals()
