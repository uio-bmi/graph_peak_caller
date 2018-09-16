import numpy as np
from .indel_scores import get_end_values_func, get_start_values_func


def max_path_func(pileup, graph, variant_maps):
    not_variant = (variant_maps.snps == 0) & (variant_maps.insertions == 0)
    get_end_value = get_end_values_func(pileup, graph)
    get_start_value = get_start_values_func(pileup, graph)

    def classify_node(prev_node, node):
        if variant_maps.deletions.from_ids[prev_node] & variant_maps.deletions.to_ids[node]:
            return "DEL"
        if variant_maps.insertions[node-graph.min_node]:
            return "INS"
        if variant_maps.snps[node-graph.min_node]:
            return "SNP"
        return "REF"

    def get_max_path(node_ids):
        node_set = set(node_ids)
        node_ids = np.asanyarray(node_ids)
        start_values = dict(zip(node_ids, get_start_value(node_ids-graph.min_node)))
        end_values = dict(zip(node_ids, get_end_value(node_ids-graph.min_node)))

        node_ids = node_ids-graph.min_node
        non_variants = node_ids[not_variant[node_ids]]
        start = np.min(non_variants) if len(non_variants) else np.min(node_ids)
        cur_node = start+graph.min_node
        path = []
        while True:
            path.append(cur_node)
            next_nodes = [node for node in graph.adj_list[int(cur_node)]
                          if node in node_set]
            if not next_nodes:
                break
            var_dict = {"REF": [], "INS": [], "DEL": [], "SNP": []}
            for node in next_nodes:
                var_dict[classify_node(cur_node, node)].append(node)
            best_vars = {key: max((start_values[node], node) for node in nodes) if len(nodes) else (-1, -1)
                         for key, nodes in var_dict.items()}
            current_best = best_vars["REF"]
            if best_vars["SNP"][0] > best_vars["REF"][0]:
                current_best = best_vars["SNP"]
            # current_best = max(best_vars["REF"], best_vars["SNP"])
            if best_vars["INS"][0] >= current_best[0]/2:
                current_best = best_vars["INS"]
            cur_value = end_values[cur_node]
            if current_best[0] <= cur_value/2:
                if current_best[0] <= best_vars["DEL"][0]/2:
                    current_best = best_vars["DEL"]
            cur_node = current_best[1]
        return path
    return get_max_path
