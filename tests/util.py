from graph_peak_caller.sparsediffs import SparseValues, SparseDiffs
from graph_peak_caller.sample.sparsegraphpileup import ReadsAdderWDirect,\
    SparseGraphPileup


def convert_old_sparse(old_sparse):
    graph = old_sparse.graph
    indices = [0]
    values = [0]
    for key in sorted(old_sparse.data.keys()):
        vi = old_sparse.data[key]
        node_idx = graph.node_indexes[key-graph.min_node]
        indices.append(node_idx)
        values.append(vi.start_value)
        indices.extend(vi.indexes+node_idx)
        values.extend(vi.values)

    sv = SparseValues(indices, values, sanitize=True)
    sv.track_size = graph.node_indexes[-1]
    return sv


def from_intervals(graph, intervals):
    pileup = SparseGraphPileup(graph)
    reads_adder = ReadsAdderWDirect(graph, pileup)
    reads_adder.add_reads(intervals)
    pileup.ends += reads_adder.pos_read_ends
    pileup.starts += reads_adder.neg_read_ends
    return SparseDiffs.from_pileup(pileup, graph.node_indexes)


def values_from_intervals(graph, intervals):
    return from_intervals(graph, intervals).get_sparse_values()
