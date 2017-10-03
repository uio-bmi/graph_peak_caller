import offsetbasedgraph as obg

Graph = obg.GraphWithReversals


def get_small_cyclic_graph():
    nodes = {1: obg.Block(100)}
    adj_list = {1: [1]}
    return Graph(nodes, adj_list)


def get_large_cyclic_graph():
    nodes = {1: obg.Block(100), 2: obg.Block(20)}
    adj_list = {1: [2], 2: [1]}
    return Graph(nodes, adj_list)


def get_padded_cyclic_graph():
    nodes = {1: obg.Block(100), 2: obg.Block(100), 3: obg.Block(100), 4: obg.Block(100)}
    adj_list = {1: [2], 2: [3, 4, 2]}
    return Graph(nodes, adj_list)
