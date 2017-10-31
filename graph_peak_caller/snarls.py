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


def build_snarlgraphs_from_vg_snarls(vg_snarls_file_name):
    snarls = Snarls.from_vg_snarls_file(vg_snarls_file_name)

    for snarl in snarls:
        if hasattr(snarl, "parent"):
            print(snarl)


