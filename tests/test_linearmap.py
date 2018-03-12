import pytest
import numpy as np
import offsetbasedgraph as obg

from graph_peak_caller.control.linearmap import LinearMap
from graph_peak_caller.control.linearintervals import\
    LinearIntervalCollection


@pytest.fixture
def hierarchical_graph():
    nodes = {i: obg.Block(i-90) for i in range(100, 106)}
    edges = {100: [101, 102],
             101: [105],
             102: [103, 104],
             103: [105],
             104: [105]}
    return obg.Graph(nodes, edges)


@pytest.fixture
def hierarchical_map():
    true_starts = [0, 10, 10, 22, 22, 36]
    true_ends = [10, 36, 22, 36, 36, 51]
    return LinearMap(true_starts, true_ends, hierarchical_graph())


def test_from_graph(hierarchical_graph, hierarchical_map):
    linear_map = LinearMap.from_graph(hierarchical_graph)
    assert linear_map == hierarchical_map


def test_get_scale_and_offset(hierarchical_map):
    node_ids = [100, 101, 102, 103, 104, 105]
    scales = [1, 26/11, 1, 14/13, 1, 1]
    offsets = [0, 10, 10, 22, 22, 36]
    for node_id, scale, offset in zip(node_ids, scales, offsets):
        assert hierarchical_map.get_scale_and_offset(node_id) == (scale, offset)


def test_map_interval_collection(hierarchical_map):
    interval_rps = [[100, 101], [100, 102],
                    [101, 105], [102, 103],
                    [102, 104], [103, 105],
                    [104, 105]]
    intervals = [obg.Interval(5, 5, rps) for rps in
                 interval_rps]
    linear_intervals = hierarchical_map.map_interval_collection(intervals)
    true_starts = [5, 5,
                   10+5*26/11, 10+5,
                   10+5, 22+5*14/13,
                   22+5]

    true_ends = [10+5*26/11, 10+5,
                 36+5, 22+5*14/13,
                 22+5, 36+5,
                 36+5]
    true_linear = LinearIntervalCollection(true_starts, true_ends)
    assert linear_intervals == true_linear
