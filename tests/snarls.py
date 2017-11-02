from graph_peak_caller.snarls import SnarlGraphBuilder
from offsetbasedgraph import Graph
from offsetbasedgraph.graphtraverser import GraphTravserserBetweenNodes
import pyvg
import sys

#vg_graph = pyvg.Graph.create_from_file("haplo1kg50-mhc.json")
#ob_graph = vg_graph.get_offset_based_graph()
#ob_graph.to_file("haplo1kg50-mhc.obg")
graph = Graph.from_file("haplo1kg50-mhc.obg")

print(graph.adj_list[153])
print(graph.adj_list[485900])
print(graph.reverse_adj_list[-485900])

#assert 598821 in graph.blocks

#traverser = GraphTravserserBetweenNodes(graph)
#subgraph = traverser.get_snarl_subgraph(153, 485900, include_start_and_end=True, print_debug=True)
#sys.exit()

#subgraph = traverser.get_greedy_subgraph_between_nodes(26577, 26579)
#print(subgraph)
#subgraph = traverser.get_greedy_subgraph_between_nodes(26571, 492107, print_debug=True)
#print(subgraph)
#sys.exit()

#subgraph = traverser.get_greedy_subgraph_between_nodes(485840, 598821)
#subgraph = traverser.get_greedy_subgraph_between_nodes(485838, 485850)
#print(subgraph)
#sys.exit()

print("N blocks: %d" % len(graph.blocks))
"""
print(graph.adj_list[598826])

print(graph.adj_list[485824])
print(graph.adj_list[485850])
print(graph.adj_list[485851])
print(graph.adj_list[485852])
print(graph.adj_list[485853])
print(graph.adj_list[485854])
"""

"""
print(graph.adj_list[126022])
print(graph.adj_list[-126022])
print(graph.adj_list[126021])
print(graph.adj_list[-126021])
print(graph.adj_list[126019])
print(graph.adj_list[-126019])
print(graph.adj_list[126020])
"""

builder = SnarlGraphBuilder.from_vg_snarls(graph, "haplo1kg50-mhc.snarls")
snarlgraph = builder.build_snarl_graphs()

#print(len(snarlgraph.blocks))

#print(len(snarlgraph.blocks[873563].blocks.keys()))


print("N blocks: %d" % len(snarlgraph.blocks))

from graph_peak_caller.snarls import SnarlGraph
from offsetbasedgraph import Block

counter = 0
counter_blocks = 0

def count_snarls(snarlgraph):
    global counter
    for block in snarlgraph.blocks.values():
        if isinstance(block, SnarlGraph):
            counter += 1
            count_snarls(block)

def count_blocks(snarlgraph):
    global counter_blocks
    for block in snarlgraph.blocks.values():
        if isinstance(block, Block):
            counter_blocks += 1
        else:
            count_blocks(block)


count_snarls(snarlgraph)
print("N snarls in snarlgraph: %d" % counter)

count_blocks(snarlgraph)
print("N blocks in snarlgraph: %d" % counter_blocks)


#for snarl in builder.snarls:
#    print(snarl)
