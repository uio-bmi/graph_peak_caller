from scipy.sparse import csr_matrix
from itertools import chain
from .fimowrapper import *
from .peakcollection import *
import numpy as np
from scipy.sparse.csgraph import shortest_path


class DiffExpression:
    def __init__(self, peak_id, main_path, var_path, main_count, var_count):
        self.peak_id = peak_id
        self.main_path = main_path
        self.var_path = var_path
        self.main_count = main_count
        self.var_count = var_count

    def __repr__(self):
        return "%s (%s)\n%s (%s)" % (self.main_path, self.main_count,
                                     self.var_path, self.var_count)


class MotifLocation:
    def __init__(self, peak, start, end):
        self._id = peak.unique_id
        self._peak = peak
        self._start = start
        self._end = end
        self.location = self._peak.get_subinterval(start, end)

    @classmethod
    def from_fimo_and_peaks(cls, fimo_entry, peaks):
        mpeaks = [p for p in peaks if p.unique_id == fimo_entry.peak_id]
        assert len(mpeaks) == 1
        return cls(mpeaks[0], fimo_entry._start-1, fimo_entry._end-1)


def main(ff, pc, subgraphs, node_ids, graph):
    peaks = list(pc.intervals)
    for i, p in enumerate(peaks):
        p.graph = graph
        p.unique_id = "peak%s" % i
    motif_locations = [MotifLocation.from_fimo_and_peaks(entry, peaks)
                       for entry in chain.from_iterable(
                               ff._entry_dict.values())]
    diffs = [get_differential(motif_location, subgraphs, node_ids)
             for motif_location in motif_locations]
    diffs = [diff for diff in diffs if diff is not None]
    for diff in diffs:
        print(diff)
    return diffs


def is_differential(motif_location, sub_graphs, node_ids_list):
    peak_id = motif_location._id
    all_rps = motif_location._peak.region_paths
    motif_rps = motif_location.location.region_paths
    graph = sub_graphs[peak_id][()]*-1
    node_ids = list(node_ids_list[peak_id])
    graph_idxs = [node_ids.index(rp) for rp in all_rps]
    motif_idxs = [node_ids.index(rp) for rp in motif_rps]
    node_scores = np.asarray(np.max(graph, axis=1).todense())
    node_scores = node_scores.flatten()
    new_data = node_scores[graph.indices]
    new_graph = csr_matrix((new_data, graph.indices, graph.indptr))
    are_differential = [i for i, row in enumerate(new_graph)
                        if i in motif_idxs and np.flatnonzero(row.data).size > 1]
    for idx in are_differential:
        in_motif = [i in motif_idxs for i in new_graph[idx].indices]
        if any(in_motif) and not all(in_motif):
            return True
    return False


def backtrace(predecessors, start, end):
    row = predecessors[start]
    cur = end
    path = [cur]
    while cur != start:
        cur = row[cur]
        path.append(cur)
    return path[::-1]


def find_alt(dists, motif_nodes, valid_nodes):
    start_node = motif_nodes[0]
    end_node = motif_nodes[-1]
    mask = valid_nodes.copy()
    mask[motif_nodes] = False
    search_nodes = np.flatnonzero(mask)
    if not len(search_nodes):
        return None, None
    start_dists = dists[start_node, search_nodes]
    end_dists = dists[search_nodes, end_node]
    tot_dists = start_dists + end_dists
    if not np.min(tot_dists) < 0:
        return None, None
    return np.flatnonzero(mask)[np.argmin(tot_dists)], tot_dists[np.argmin(tot_dists)]


def get_differential(motif_location, sub_graphs, node_ids_list):
    peak_id = motif_location._id
    motif_rps = motif_location.location.region_paths
    graph = sub_graphs[peak_id][()]
    node_ids = list(node_ids_list[peak_id])
    motif_idxs = [node_ids.index(rp) for rp in motif_rps]
    node_scores = np.asarray(np.min(graph, axis=1).todense())
    node_scores = node_scores.flatten()
    valid_nodes = node_scores != 0
    dists, predecessors = shortest_path(graph, return_predecessors=True,
                                        method="BF")
    alt, altD = find_alt(dists, motif_idxs, valid_nodes)
    if alt is None:
        return None
    mainD = dists[motif_idxs[0], motif_idxs[-1]]
    dist_diff = altD-mainD
    alt_node_score = abs(np.min(graph[alt]))
    start_path = backtrace(predecessors, motif_idxs[0], alt)
    end_path = backtrace(predecessors, alt, motif_idxs[-1])
    path = start_path[:-1] + end_path
    node_idxs = node_ids_list[peak_id][path]
    alt_interval = obg.DirectedInterval(
        motif_location.location.start_position.offset,
        motif_location.location.end_position.offset,
        list(node_idxs))
    obg_graph = motif_location.location.graph
    alt_interval.graph = obg_graph
    length_diff = motif_location.location.length()-alt_interval.length()
    alt_length = obg_graph.node_size(node_ids[alt])
    base_score = alt_node_score/alt_length
    main_score = base_score + dist_diff
    return DiffExpression(peak_id, motif_location.location,
                          alt_interval,
                          main_score, base_score)
