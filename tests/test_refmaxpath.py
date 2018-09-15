import offsetbasedgraph as obg
from collections import defaultdict
from graph_peak_caller.postprocess.reference_based_max_path\
    import max_path_func
from offsetbasedgraph.vcfmap import VariantMap, DELMap
from graph_peak_caller.sparsediffs import SparseValues
from graph_peak_caller.peakcollection import Peak
from graph_peak_caller.postprocess.maxpaths import SparseMaxPaths
import numpy as np
import pytest


@pytest.fixture
def graph():
    nodes = {i: obg.Block(10) for i in range(10, 18)}
    edges = {10: [11], 11: [12, 13], 12: [14], 13: [14],
             14: [15, 16, 17], 15: [16], 16: [17]}
    return obg.Graph(nodes, edges)


@pytest.fixture
def variant_maps():
    snps, insertions = np.zeros((2, 8))
    snps[3] = 1
    insertions[5] = 1
    from_ids, to_ids = (defaultdict(set), defaultdict(set))
    from_ids[14].add(1)
    to_ids[17].add(1)
    deletions = DELMap(from_ids, to_ids)
    return VariantMap(snps=snps,
                      insertions=insertions,
                      deletions=deletions,
                      trails=np.zeros_like(snps))


@pytest.fixture
def pileup():
    sv = SparseValues([0, 15, 30, 45, 50, 60, 70],
                      [0, 1,   2,  6,  2, 3, 6])
    sv.track_size = 80
    return sv


@pytest.fixture
def nonvar_pileup():
    sv = SparseValues([0, 15, 45, 50, 60, 70],
                      [0, 1,   4,  1,  3, 5])
    sv.track_size = 80
    return sv


@pytest.fixture
def get_max_path():
    return max_path_func(pileup(), graph(), variant_maps())


@pytest.fixture
def get_max_path_nonvar():
    return max_path_func(nonvar_pileup(), graph(), variant_maps())


@pytest.fixture
def areas():
    sv = SparseValues([0, 15, 42, 48,  72],
                      [0, 1, 0, 1, 0])
    sv.track_size = 80
    return sv


@pytest.fixture
def offset_areas():
    sv = SparseValues([0, 10, 11, 20, 29, 31],
                      [1, 0, 1, 0, 1, 0])
    sv.track_size = 80
    return sv


@pytest.fixture
def sparse_max_paths():
    return SparseMaxPaths(areas(), graph(), pileup(), variant_maps())


@pytest.fixture
def sparse_max_paths_nonvar():
    return SparseMaxPaths(areas(), graph(), nonvar_pileup(), variant_maps())


@pytest.fixture
def sparse_max_paths_offset():
    return SparseMaxPaths(offset_areas(), graph(), nonvar_pileup(), variant_maps())


def test_snp(get_max_path):
    max_path = get_max_path([11, 12, 13, 14])
    assert max_path == [11, 13, 14]


def test_ins(get_max_path):
    max_path = get_max_path([14, 15, 16])
    assert max_path == [14, 15, 16]


def test_del(get_max_path):
    max_path = get_max_path([14, 15, 16, 17])
    assert max_path == [14, 17]


def test_no_snp(get_max_path_nonvar):
    max_path = get_max_path_nonvar([11, 12, 13, 14])
    assert max_path == [11, 12, 14]


def test_no_ins(get_max_path_nonvar):
    max_path = get_max_path_nonvar([14, 15, 16])
    assert max_path == [14, 16]


def test_no_del(get_max_path_nonvar):
    max_path = get_max_path_nonvar([14, 15, 16, 17])
    assert max_path == [14, 16, 17]


def test_maxpath(sparse_max_paths):
    max_paths = sparse_max_paths.run()[0]
    print(max_paths)
    true_max_paths = [Peak(5, 2, [11, 13, 14]),
                      Peak(8, 2, [14, 17])]
    max_paths.sort(key=lambda peak: peak.region_paths[0])
    assert max_paths == true_max_paths


def test_maxpath_novar(sparse_max_paths_nonvar):
    max_paths = sparse_max_paths_nonvar.run()[0]
    true_max_paths = [Peak(5, 2, [11, 12, 14]),
                      Peak(8, 2, [14, 16, 17])]
    max_paths.sort(key=lambda peak: peak.region_paths[0])
    assert max_paths == true_max_paths


def test_start_offset(sparse_max_paths_offset):
    offsets =  sparse_max_paths_offset.get_start_offsets([10, 11, 12])
    assert list(offsets) == [0, 1, 9]


def test_end_offset(sparse_max_paths_offset):
    sparse_max_paths_offset._sparse_values.values = 1-sparse_max_paths_offset._sparse_values.values
    offsets =  sparse_max_paths_offset.get_end_offsets([11, 12])
    assert list(offsets) == [1, 9]


def test_old_maxpath(sparse_max_paths):
    sparse_max_paths._variant_maps = None
    max_paths = sparse_max_paths.run()[0]
    true_max_paths = [Peak(5, 2, [11, 13, 14]),
                      Peak(8, 2, [14, 15, 16, 17])]
    max_paths.sort(key=lambda peak: peak.region_paths[0])
    assert max_paths == true_max_paths


#    sv = SparseValues([0, 10, 11, 20, 29, 31],
