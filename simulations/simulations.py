
from graphsimulator import GraphSimulator
from pileupsimulator import PileupSimulator
from graph_peak_caller.callpeaks import CallPeaks, ExperimentInfo
from graph_peak_caller.snarls import SnarlGraphBuilder
from graph_peak_caller.linearsnarls import LinearSnarlMap
from offsetbasedgraph import IntervalCollection
from graph_peak_caller.peakcollection import PeakCollection
import logging
logging.basicConfig(level=logging.INFO)


class SimulatedPeakCalling():
    def __init__(self, n_paths, n_basepairs_length, n_snps, n_peaks, with_control=False):
        self.n_paths = n_paths
        self.n_basepairs_length = n_basepairs_length
        self.n_snps = n_snps
        self.n_peaks = n_peaks
        self.snarls = None
        self.graph = None
        self.sample_pileup = None
        self.control_pileup = None
        self.correct_peaks = None
        self.n_sample_reads = None
        self.n_control_reads = None
        self.control_reads = None
        self.with_control = with_control

        self.set_up()

    def set_up(self):
        graph_simulator = GraphSimulator(n_paths=self.n_paths,
                                n_basepairs_length=self.n_basepairs_length,
                                n_snps=self.n_snps)
        simulated_graph = graph_simulator.get_simulated_graph()
        self.snarls = simulated_graph.snarls

        pileup_simulator = PileupSimulator(
                                simulated_graph=simulated_graph,
                                n_peaks = self.n_peaks,
                                with_control=self.with_control)

        self.correct_peaks = pileup_simulator.get_correct_peak_positions_on_graph()
        sample_pileup, control_pileup = pileup_simulator.get_simulated_pileups()
        self.sample_pileup = sample_pileup
        self.control_pileup = control_pileup
        self.graph = simulated_graph.graph.to_graph_with_reversals()
        self.n_sample_reads = pileup_simulator.n_sample_reads
        self.n_control_reads = pileup_simulator.n_control_reads
        self.control_reads = pileup_simulator._control_reads

        print("Graph")
        print(self.graph)

    def call_peaks(self):

        genome_size = sum(block.length() for block in self.graph.blocks.values())
        experiment_info = ExperimentInfo(genome_size, 50, 20)
        experiment_info.n_sample_reads = self.n_sample_reads
        experiment_info.n_control_reads = self.n_control_reads


        snarlbuilder = SnarlGraphBuilder(self.graph.copy(), self.snarls,
                                         id_counter=self.graph.max_block_id() + 1)
        snarlgraph = snarlbuilder.build_snarl_graphs()
        linear_map = LinearSnarlMap(snarlgraph, self.graph)
        linear_map.to_file("simulated_snarl_map.tmp")
        caller = CallPeaks(self.graph,
                           sample_intervals="dummy",
                           control_intervals=IntervalCollection(self.control_reads),
                           experiment_info=experiment_info,
                           has_control=self.with_control,
                           linear_map="simulated_snarl_map.tmp")

        caller._sample_pileup = self.sample_pileup
        caller.create_control()
        caller.scale_tracks()
        caller.get_score()
        caller.call_peaks("simulated.peaks")

    def compare_with_correct_peaks(self):
        correct_peaks = PeakCollection(self.correct_peaks)
        for peak in correct_peaks:
            print(peak)
        found_peaks = PeakCollection.create_list_from_file("max_paths", graph=self.graph)

        for i in found_peaks:
            print(i)
        matched = correct_peaks.get_identical_intervals(found_peaks)

        print("%d correct peaks identically found, %3.f %% " % (len(matched), 100 * len(matched) / len(correct_peaks.intervals)))
        for i in correct_peaks:
            print(i)


if __name__ == "__main__":
    """
    simulator = GraphSimulator(2, 1000, 25)
    simulated_graph = simulator.get_simulated_graph()
    print(simulated_graph.graph)
    import sys
    sys.exit()
    """
    caller = SimulatedPeakCalling(
        n_paths=2,
        n_basepairs_length=10000,
        n_snps = 10,
        n_peaks = 5,

        with_control=False
    )

    #print(caller.graph)
    caller.call_peaks()
    caller.compare_with_correct_peaks()
