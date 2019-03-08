import pytest
import numpy as np
from offsetbasedgraph.vcfgraph import VCFGraph, Interval, AdjList, SNPs
from graph_peak_caller.vcfpileup import ForwardPileup, ReversePileup, MainPileup
from graph_peak_caller.linearrepr import LinearRepresentation
from graph_peak_caller.sparsediffs import SparseDiffs


@pytest.fixture
def graph():
    node_lens = [10, 10, 10]
    adj_list = {0: [1, 2], 1: [2]}
    return VCFGraph(node_lens, AdjList.from_dict(adj_list, 3), SNPs())


@pytest.fixture
def pos_intervals():
    return [Interval(5, 7, [0], 1),
            Interval(3, 5, [1], 1),
            Interval(1, 3, [2], 1)]


@pytest.fixture
def neg_intervals():
    return [Interval(5, 7, [0], -1),
            Interval(3, 5, [1], -1),
            Interval(1, 3, [2], -1)]


def test_forward_extend(graph, pos_intervals):
    node_counter = np.zeros((graph._node_lens.size, 2), dtype="int")
    pileup = ForwardPileup(graph, 5, node_counter)
    for interval in pos_intervals:
        pileup.add_interval(interval)
    assert np.all(node_counter == np.zeros_like(node_counter))
    pileup.extend_intervals()
    true_node_counts = np.array([[0, 1],
                                 [1, 0],
                                 [1, 0]])
    assert np.all(node_counter == true_node_counts)
    true_starts = [[5], [3], [1]]
    true_ends = [[], [2, 10], [2, 8]]
    for res, true in zip(pileup.starts, true_starts):
        assert list(res) == true

    for res, true in zip(pileup.ends, true_ends):
        assert list(res) == true


def test_long_extend(graph, pos_intervals):
    node_counter = np.zeros((graph._node_lens.size, 2), dtype="int")
    pileup = ForwardPileup(graph, 10, node_counter)
    for interval in pos_intervals:
        pileup.add_interval(interval)
    assert np.all(node_counter == np.zeros_like(node_counter))
    pileup.extend_intervals()
    true_node_counts = np.array([[0, 1],
                                 [1, 1],
                                 [2, 1]])
    assert np.all(node_counter == true_node_counts)
    true_starts = [[5], [3], [1]]
    true_ends = [[], [7], [5, 7]]
    for res, true in zip(pileup.starts, true_starts):
        assert list(sorted(res)) == true

    for res, true in zip(pileup.ends, true_ends):
        assert list(sorted(res)) == true


def test_reverse_extend(graph, neg_intervals):
    node_counter = np.zeros((graph._node_lens.size, 2), dtype="int")
    pileup = ReversePileup(graph, 5, node_counter)
    for interval in neg_intervals:
        pileup.add_interval(interval)
    assert np.all(node_counter == np.zeros_like(node_counter))
    pileup.extend_intervals()
    true_node_counts = np.array([[0, 2],
                                 [1, 1],
                                 [1, 0]])
    assert np.all(node_counter == true_node_counts)
    true_starts = [[0, 6, 8], [6], []]
    true_ends = [[7], [5], [3]]
    for res, true in zip(pileup.starts, true_starts):
        assert list(sorted(res)) == true

    for res, true in zip(pileup.ends, true_ends):
        assert list(sorted(res)) == true


def test_long_reverse_extend(graph, neg_intervals):
    node_counter = np.zeros((graph._node_lens.size, 2), dtype="int")
    pileup = ReversePileup(graph, 10, node_counter)
    for interval in neg_intervals:
        pileup.add_interval(interval)
    assert np.all(node_counter == np.zeros_like(node_counter))
    pileup.extend_intervals()
    true_node_counts = np.array([[1, 2],
                                 [1, 1],
                                 [1, 0]])
    assert np.all(node_counter == true_node_counts)
    true_starts = [[1, 3], [1], []]
    true_ends = [[7], [5], [3]]
    for res, true in zip(pileup.starts, true_starts):
        assert list(sorted(res)) == true

    for res, true in zip(pileup.ends, true_ends):
        assert list(sorted(res)) == true


def test_pileup(graph, pos_intervals, neg_intervals):
    lin_rep = LinearRepresentation(graph)
    pileup = MainPileup(graph, lin_rep, 5)
    sd = pileup.build(pos_intervals+neg_intervals)
    true_indexes = [0, 5, 6,  7, 8, 10, 12, 13, 15, 16, 20, 21, 22, 23, 28]
    true_values =  [1, 1, 1, -1, 1, -1, -1,  1, -1,  1,  0,  1, -1, -1, -1]
    true_sd = SparseDiffs(true_indexes, true_values)
    assert sd.get_sparse_values() == true_sd.get_sparse_values()
