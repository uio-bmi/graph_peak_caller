import pytest
import offsetbasedgraph as obg
from graph_peak_caller.linear_filter import LinearFilter


def graph():
    nodes = [9, 10, 11, 12, 13, 14]
    nodes = {i: obg.Block(10) for i in nodes}
    edges = {
        9: [10],
        10: [11, 12],
        11: [13],
        12: [13],
        13: [14]}
    return obg.Graph(nodes, edges)


@pytest.fixture
def indexed_interval():
    start = obg.Position(10, 5)
    end = obg.Position(13, 5)
    return obg.IndexedInterval(
        start, end, [10, 11, 13],
        graph=graph())


def test_get_start_positions(indexed_interval):
    positions = [obg.Position(10, 7),
                 obg.Position(11, 5),
                 obg.Position(-11, 5),
                 obg.Position(13, 2)]
    linear_filter = LinearFilter(positions, indexed_interval)
    start_positions = linear_filter.find_start_positions()
    assert start_positions == {"+": [2, 10, 17],
                               "-": [10]}
