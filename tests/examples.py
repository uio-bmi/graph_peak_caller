import numpy as np
import json
import offsetbasedgraph
from graph_peak_caller.vg import *
from graph_peak_caller.pileup import Pileup

position_jsons = [json.loads(pos_str) for pos_str in
                  ('{"offset": 10, "node_id": 0}',
                   '{"offset": 0, "node_id": 1, "is_reverse": true}')]

positions = [Position(0, 10, False), Position(1, 0, True)]

edit_jsons = [json.loads(edit_str) for edit_str in
              ['{"to_length": 4, "from_length": 4}',
               '{"to_length": 1, "from_length": 1, "sequence": "N"}']]

edits = [Edit(4, 4, None), Edit(1, 1, "N")]

mapping_jsons = [
    {"position": position_jsons[0], "edit": edit_jsons},
    {"position": position_jsons[1], "edit": edit_jsons}]

mappings = [Mapping(positions[0], edits),
            Mapping(positions[1], edits)]

path_jsons = [
    {"name": "path1",
     "mapping": mapping_jsons[:1]}]

paths = [Path("path1", mappings[:1])]
intervals = [offsetbasedgraph.Interval(10, 15, [0])]

nodes = [Node("_", 0, 10),
         Node("_", 1, 20),
         Node("_", 2, 30),
         Node("_", 3, 40)]

edges = [Edge(0, 1),
         Edge(0, 2),
         Edge(1, 3),
         Edge(2, 3)]

vg_graphs = [Graph(nodes, edges, [])]
ob_graphs = [graph.get_offset_based_graph() for graph in vg_graphs]

pileup_intervals = [offsetbasedgraph.Interval(5, 5, [0, 1]),
                    offsetbasedgraph.Interval(0, 10, [1]),
                    offsetbasedgraph.Interval(10, 15, [1])]
for interval in pileup_intervals:
    interval.graph = ob_graphs[0]

true_counts = {node.id: np.zeros(node.n_basepairs, dtype="int32")
               for node in nodes}
true_counts[0][5:] = 1
true_counts[1][:5] = 2
true_counts[1][5:15] = 1
true_pileup = Pileup(None, None)
true_pileup.count_arrays = true_counts


one_block_graph = offsetbasedgraph.Graph({1: offsetbasedgraph.Block(10)}, {})
one_block_pileup1_intervals = [offsetbasedgraph.Interval(1, 10, [1], one_block_graph),
                               offsetbasedgraph.Interval(1, 5, [1], one_block_graph)]
one_block_pileup2_intervals = [offsetbasedgraph.Interval(0, 5, [1], one_block_graph)]

pileup1_one_block = Pileup(one_block_graph, one_block_pileup1_intervals).create()
pileup2_one_block = Pileup(one_block_graph, one_block_pileup2_intervals).create()
