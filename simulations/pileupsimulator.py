from pyvg.sequences import SequenceRetriever
import offsetbasedgraph as obg
import pyvg
from graph_peak_caller.util import fasta_sequence_to_linear_path_through_graph, get_linear_paths_in_graph
import random
import numpy as np
from math import floor
from graph_peak_caller.sparsepileup import SparsePileup
from offsetbasedgraph import Interval
from graphsimulator import GraphSimulator
random.seed(1)

class PileupSimulator():

    def __init__(self, simulated_graph, n_peaks, with_control=False):

        self.simulated_graph = simulated_graph
        self.n_peaks = n_peaks
        self.peak_size = 50
        self.n_reads_at_peak = 10
        self.with_control = with_control

        self._linear_reads = []
        self._sample_reads = []
        self._control_reads = []
        self.linear_peaks = []
        self.n_sample_reads = 0
        self.n_control_reads = 0

        self._create_linear_reads()
        self._add_noise_to_linear_reads()
        self._translate_reads_to_graph()

    def get_correct_peak_positions_on_graph(self):
        peaks = []
        for peak in self.linear_peaks:
            peaks.append(self.simulated_graph.translate(peak))
        return peaks

    def _create_intervals_around_peak_position(self, node_id, offset):

        linear_interval = Interval(int(offset - self.peak_size/2),
                                   int(offset + self.peak_size/2),
                                   [node_id])
        for i in range(0, self.n_reads_at_peak):
            print("Created linear interval: %s" % linear_interval)
            #graph_interval = self.simulated_graph.translate(linear_interval)
            self._linear_reads.append(linear_interval.copy())
            self.n_sample_reads += 1

        self.linear_peaks.append(linear_interval)

    def _create_linear_reads(self):

        peak_locations = self.peak_size * np.array(
                random.sample(
                    range(1, floor((self.simulated_graph.linear_graph_length - self.peak_size)
                                   / self.peak_size), 2), self.n_peaks))

        for location in peak_locations:
            path = random.sample(range(0, self.simulated_graph.n_linear_paths), 1)
            path = 0
            path_id= path + 100
            self._create_intervals_around_peak_position(path_id, location)

    def _add_noise_to_linear_reads(self):
        pass

    def _translate_reads_to_graph(self):
        for interval in self._linear_reads:
            graph_interval = self.simulated_graph.translate(interval)
            self._sample_reads.append(graph_interval)
            if not self.with_control:
                self._control_reads.append(graph_interval)
                self.n_control_reads += 1

    def get_simulated_pileups(self):
        sample_pileup = SparsePileup.from_intervals(self.simulated_graph.graph,
                                                    self._sample_reads)
        control_pileup = SparsePileup.from_intervals(self.simulated_graph.graph,
                                                     self._control_reads)
        return sample_pileup, control_pileup

if __name__ == "__main__":

    simulator = GraphSimulator(n_paths=2,
                               n_basepairs_length=200,
                               n_snps=1)
    simulated_graph = simulator.get_simulated_graph()
    print(simulator.graph)

    pileup_simulator = PileupSimulator(
                            simulated_graph=simulated_graph,
                            n_peaks = 1,
                            with_control=False)

    linear_pileup, control_pileup = pileup_simulator.get_simulated_pileups()

    print(linear_pileup)
"""

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

        ob_graph.adj_list[6510] = [6511]
        ob_graph.reverse_adj_list[-6511] = [6510]

        print(ob_graph.adj_list[-6511])
        print(ob_graph.reverse_adj_list[-6510])
        print(ob_graph.reverse_adj_list[-6510])
        print(ob_graph.reverse_adj_list[6510])
        print(ob_graph.adj_list[6510])
        #return
        print(ob_graph.reverse_adj_list[-5703])


        self.get_linear_paths()

    def get_linear_paths(self):
        for file_name, start in [(self.path2_sequence_file_name, self.path2_start_node),
                                 (self.path1_sequence_file_name, self.path1_start_node)]:
            print("Finding path for %s" % file_name)
            sequence = open(file_name).read()
            #print(sequence[48385-500:48385 + 500])
            #print(sequence[46094:46094 + 4000])
            #return

            linear_path = fasta_sequence_to_linear_path_through_graph(file_name, self.sequence_retriever, self.ob_graph, start)
            self.linear_paths.append(linear_path)
            print("Found linear path. Length: %d" % linear_path.length())



vg_graph = vg_graph = pyvg.Graph.create_from_file("../tests/cactus-mhc.json")
ob_graph = obg.GraphWithReversals.from_file("../tests/cactus-mhc.obg")
get_linear_paths_in_graph(ob_graph, vg_graph, "linear_paths.intervalcollection")

"""



"""
sequence_retriever = SequenceRetriever.from_vg_graph("../tests/cactus-mhc.vg")
ob_graph = obg.GraphWithReversals.from_file("../tests/cactus-mhc.obg")

vg_graph = vg_graph = pyvg.Graph.create_from_file("../tests/cactus-mhc.json")
#ob_graph = vg_graph.get_offset_based_graph()

#from pyvg.vg import ProtoGraph
#vg_graph = ProtoGraph.from_vg_graph_file("../tests/cactus-mhc.vg", only_read_nodes=False)


simulator = PileupSimulator(ob_graph, sequence_retriever,
                            "mhc_cleaned2.fa", 225518,
                            "GI568335992_cleaned.fa", 5619,
                            100
                            )

"""
"""
simulator = PileupSimulator(ob_graph, sequence_retriever,
                            "mhc_cleaned2.fa", 225518,
                            "sample.fasta", 6477,
                            100
                            )
"""



#GI568335986_cleaned


