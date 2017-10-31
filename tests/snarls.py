from graph_peak_caller.snarls import SnarlGraphBuilder
from offsetbasedgraph import Graph

graph = Graph.from_file("cactus-mhc.obg")
builder = SnarlGraphBuilder.from_vg_snarls(graph, "snarls.pb")
builder.build_snarl_graphs()

#for snarl in builder.snarls:
#    print(snarl)
