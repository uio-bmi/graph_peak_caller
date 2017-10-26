from pyvg.sequences import SequenceRetriever
import offsetbasedgraph as obg
import pyvg
from graph_peak_caller.util import fasta_sequence_to_linear_path_through_graph

class PileupSimulator(object):

    def __init__(self, ob_graph, sequence_retriever,
                        path1_sequence_file_name,
                        path1_start_node,
                        path2_sequence_file_name,
                        path2_start_node,
                        n_peaks):
        self.ob_graph = ob_graph
        self.sequence_retriever = sequence_retriever
        self.path1_sequence_file_name = path1_sequence_file_name
        self.path2_sequence_file_name = path2_sequence_file_name
        self.path1_start_node = path1_start_node
        self.path2_start_node = path2_start_node
        self.n_peaks = n_peaks
        self.linear_paths = []

        self.get_linear_paths()

    def get_linear_paths(self):
        for file_name, start in [(self.path2_sequence_file_name, self.path2_start_node),
                                 (self.path1_sequence_file_name, self.path1_start_node)]:
            print("Finding path for %s" % file_name)
            linear_path = fasta_sequence_to_linear_path_through_graph(file_name, self.sequence_retriever, self.ob_graph, start)
            self.linear_paths.append(linear_path)
            print("Found linear path. Length: %d" % linear_path.length())





if __name__ == "__main__":
    sequence_retriever = SequenceRetriever.from_vg_graph("../tests/cactus-mhc.vg")
    ob_graph = obg.GraphWithReversals.from_file("../tests/cactus-mhc.obg")
    #vg_graph = vg_graph = pyvg.Graph.create_from_file("../testscactus-mhc.json")
    simulator = PileupSimulator(ob_graph, sequence_retriever,
                                "mhc_cleaned2.fa", 225518,
                                "GI568335986_cleaned.fa", 5619,
                                100
                                )


