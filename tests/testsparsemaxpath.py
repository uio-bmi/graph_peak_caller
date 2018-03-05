import offsetbasedgraph as obg
from graph_peak_caller.postprocess.maxpaths import SparseMaxPaths
from graph_peak_caller.sparsediffs import SparseValues
from graph_peak_caller.peakcollection import Peak
import numpy as np


nodes = {i+1: obg.Block(10) for i in range(10)}
edges = {i: [i+1] for i in range(1, 10)}
edges[1] = [2, 3]
edges[2] = [4]
edges[3] = [4]

graph = obg.GraphWithReversals(nodes, edges)


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


def offsetgraph():
    sizes = np.ones(11)*10
    sizes[0] = 0
    blocks = obg.BlockArray(sizes)
    blocks.node_id_offset = 100
    # nodes = {i+1: obg.Block(10) for i in range(100, 110)}
    edges = {i: [i+1] for i in range(101, 110)}
    edges[101] = [102, 103]
    edges[102] = [104]
    edges[103] = [104]
    return obg.GraphWithReversals(blocks, edges)


def test_single_peak():
    indices = np.array([0, 3, 32], dtype="int")
    values = np.array([False, True, False], dtype="bool")
    pileup = SparseValues(indices, values)
    indices = np.arange(0, 40, 2)
    values = (np.arange(0, 40, 2)+1) % 5
    score_pileup = SparseValues(indices, values)
    score_pileup.track_size = 100
    graph.node_indexes = np.arange(0, 110, 10)
    max_paths = SparseMaxPaths(pileup, graph, score_pileup)
    print(max_paths.run())


def test_simple_peak():
    pileup = SparseValues([0, 5, 35], [False, True, False])
    score_pileup = SparseValues([0, 5, 10, 20, 30, 35],
                                [0, 1, 2, 3, 4, 0])
    score_pileup.track_size = 100
    pileup.track_size = 100
    max_paths = SparseMaxPaths(pileup, graph, score_pileup).run()
    assert max_paths == [Peak(5, 5, [1, 3, 4], graph=graph)]


def test_offset_peak():
    graph = offsetgraph()
    pileup = SparseValues([0, 5, 35], [False, True, False])
    score_pileup = SparseValues([0, 5, 10, 20, 30, 35],
                                [0, 1, 2, 3, 4, 0])
    score_pileup.track_size = 100
    pileup.track_size = 100
    max_paths = SparseMaxPaths(pileup, graph, score_pileup).run()
    print(max_paths)
    assert max_paths == [Peak(5, 5, [101, 103, 104], graph=graph)]


def test_trailing_zeros():
    graph = offsetgraph()
    pileup = SparseValues([0, 5, 35], [False, True, False])
    score_pileup = SparseValues([0, 15, 30],
                                [0, 10, 0])
    score_pileup.track_size = 100
    pileup.track_size = 100
    max_paths = SparseMaxPaths(pileup, graph, score_pileup).run()
    print(max_paths)
    assert max_paths == [Peak(5, 5, [101, 103, 104], graph=graph)]


def test_offset_end_peak():
    graph = offsetgraph()
    pileup = SparseValues([0, 85], [False, True])
    score_pileup = SparseValues([0, 85],
                                [0, 10])
    score_pileup.track_size = 100
    pileup.track_size = 100
    max_paths = SparseMaxPaths(pileup, graph, score_pileup).run()
    print(max_paths)
    assert max_paths == [Peak(5, 10, [109, 110], graph=graph)]


def test_multiple_peak():
    graph = complicated_offset()
    pileup = SparseValues([0, 1, 7, 8, 10, 12, 14],
                          [0, 1, 0, 1,  0,  1, 0])
    score_pileup = SparseValues([0, 4],
                                [5, 6])
    score_pileup.track_size = 20
    pileup.track_size = 20
    max_paths = SparseMaxPaths(pileup, graph, score_pileup).run()
    max_paths.sort(key=lambda x: x.region_paths[0])
    print(max_paths)
    assert max_paths == [Peak(1, 1, [101, 103, 104], graph=graph),
                         Peak(0, 2, [105, 107], graph=graph)]


if __name__ == "__main__":
    # test_simple_peak()
    # test_offset_peak()
    # test_offset_end_peak()
    # test_trailing_zeros()
    test_multiple_peak()
