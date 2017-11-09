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

        self._sample_linear_reads = []
        self._control_linear_reads = []
        self._sample_reads = []
        self._control_reads = []
        self.linear_peaks = []
        self.n_sample_reads = 0
        self.n_control_reads = 0

        self._create_linear_reads()
        self._add_noise_to_linear_reads()

        if self.with_control:
            #self.create_control_on_other_path()
            self.add_background_control_reads()
        else:
            self.copy_sample_to_control()

        self._translate_reads_to_graph()

        print("N sample reads: %d" % self.n_sample_reads)
        print("N control reads: %d" % self.n_control_reads)

    def write_reads_to_bed(self):
        f = open("sample.bed", "w")
        for peak in self._sample_linear_reads:
            line = "chr1\t%d\t%d\t.\t0\t+\n" % (peak.start_position.offset, peak.end_position.offset)
            f.writelines([line])
        f.close()
        print("Wrote to bed")

    def create_control_on_other_path(self):
        print("Creating control on other path")
        i = 0
        for read in self._sample_linear_reads:
            i += 1
            if i % 2 == 1:
                continue
            path = read.region_paths[0]
            other_path = (path - 100  + 1 % self.simulated_graph.n_linear_paths) + 100
            control_read = Interval(read.start_position.offset,
                                    read.end_position.offset,
                                    [other_path])
            self.n_control_reads += 1
            self._control_linear_reads.append(control_read)

    def _random_path_id(self):
        return 100 + random.sample(range(0, self.simulated_graph.n_linear_paths), 1)[0]

    def add_background_control_reads(self):
        #n_to_add = int(self.n_peaks * 10 * 2)
        n_to_add = int(1 * self.simulated_graph.linear_graph_length / self.peak_size)
        peak_locations = random.sample(range(self.peak_size,
                                             self.simulated_graph.linear_graph_length - self.peak_size
                                            ), n_to_add)
        for location in peak_locations:
            path = self._random_path_id()
            interval = Interval(location, location + self.peak_size, [path])
            self._control_linear_reads.append(interval)
            self.n_control_reads += 1

    def copy_sample_to_control(self):
        for read in self._sample_linear_reads:
            self._control_linear_reads.append(read)
            self.n_control_reads += 1

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
            self._sample_linear_reads.append(linear_interval.copy())
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
            path_id = path + 100
            self._create_intervals_around_peak_position(path_id, location)

    def _add_noise_to_linear_reads(self):
        pass

    def _translate_reads_to_graph(self):
        for interval in self._sample_linear_reads:
            graph_interval = self.simulated_graph.translate(interval)
            self._sample_reads.append(graph_interval)

        for interval in self._control_linear_reads:
            graph_interval = self.simulated_graph.translate(interval)
            self._control_reads.append(graph_interval)

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



