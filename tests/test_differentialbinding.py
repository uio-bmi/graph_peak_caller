import pytest
import numpy as np
from scipy.sparse import csr_matrix
import offsetbasedgraph as obg

from graph_peak_caller.analysis.differentialbinding import *
from graph_peak_caller.peakcollection import Peak


@pytest.fixture
def graph():
    nodes = {i: obg.Block(5) for i in range(1, 5)}
    edges = {1: [2, 3], 2: [4], 3: [4]}
    return obg.Graph(nodes, edges)


@pytest.fixture
def motif_location():
    peak = Peak(1, 4, [1, 3, 4], graph(), unique_id="peak1")
    return MotifLocation(peak, 2, 11)


@pytest.fixture
def subgraphs():
    subgraph = csr_matrix(([-10, -10, -20, -30, -40],
                           ([0,    0,   1,   2,   3],
                            [1,    2,   3,   3,   4])),
                          shape=[5, 5])
    return {"peak1": np.array(subgraph)}


@pytest.fixture
def node_ids():
    return {"peak1": np.array([1, 2, 3, 4], dtype="int")}


@pytest.fixture
def diff_expr():
    main_path = obg.Interval(3, 2, [1, 3, 4])
    var_path = obg.Interval(3, 2, [1, 2, 4])
    return DiffExpression("peak1", main_path, var_path,
                          14, 4)


def test_is_differential(motif_location, subgraphs, node_ids, diff_expr):
    print(subgraphs["peak1"][()])
    diff = get_differential(motif_location, subgraphs, node_ids)
    assert diff == diff_expr
