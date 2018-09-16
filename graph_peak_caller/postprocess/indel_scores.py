import numpy as np
from offsetbasedgraph.vcfmap import next_node_func, prev_node_func


def get_end_values_func(pileup, graph):
    node_indexes = graph.node_indexes

    def get_end_values(node_ids):
        node_idxs = node_indexes[node_ids+1]-1
        pileup_idxs = np.digitize(node_idxs, pileup.indices)-1
        return pileup.values[pileup_idxs]

    return get_end_values


def get_start_values_func(pileup, graph):
    node_indexes = graph.node_indexes

    def get_start_values(node_ids):
        node_idxs = node_indexes[node_ids]
        pileup_idxs = np.digitize(node_idxs, pileup.indices)-1
        return pileup.values[pileup_idxs]

    return get_start_values


def find_invalid_insertions(pileup, insertion_map, graph, reference):
    get_next_ref = next_node_func(graph, reference)
    get_prev_ref = prev_node_func(graph, reference)
    get_end_value = get_end_values_func(pileup, graph)
    get_start_value = get_start_values_func(pileup, graph)


    insertion_nodes = np.flatnonzero(insertion_map) + graph.min_node
    prev_nodes = np.array([get_prev_ref(node)-graph.min_node for node in insertion_nodes])
    next_nodes = np.array([get_next_ref(node)-graph.min_node for node in insertion_nodes])
    last_values = get_end_value(prev_nodes)
    first_values = get_start_value(next_nodes)
    avg_values = (get_start_value(insertion_nodes)+get_end_value(insertion_nodes))/2
    valid_nodes = avg_values*2 > (first_values+last_values)/2
    in_valid_map = insertion_map.copy() > 0
    print("Insertions: %s" % np.count_nonzero(insertion_map))
    in_valid_map[insertion_nodes[valid_nodes]] = 0
    print("InvalidInsertions: %s" % np.count_nonzero(in_valid_map))
    return in_valid_map

