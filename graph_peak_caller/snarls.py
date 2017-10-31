from offsetbasedgraph import Graph
from pyvg import Snarls

class SnarlGraph(Graph):

    def __init__(blocks, edges):
        pass

    @classmethod
    def create_from_start_and_end(cls, start, end, parent_graph):
        pass

    def length(self):
        pass


class SimpleSnarl():
    def __init__(self, start, end, id, parent=None):
        self.start = start
        self.end = end
        self.id = id
        self.parent = parent
        

class SnarlGraphBuilder():

    def __init__(self, graph, snarls):
        self.snarls = snarls

    @classmethod
    def from_vg_snarls(cls, vg_snarls_file_name):

        snarls = Snarls.from_vg_snarls_file(vg_snarls_file_name)

        for snarl in snarls:
            if hasattr(snarl, "parent"):
                print(snarl)


