from itertools import chain
import numpy as np
import offsetbasedgraph as obg
from pyvg.sequences import SequenceRetriever
from graph_peak_caller.callpeaks import CallPeaks, ExperimentInfo
import random
import logging
logging.basicConfig(level=logging.INFO)


class IntervalSimulator(object):
    def __init__(self, graph, interval_length):
        self._graph = graph
        self._interval_length = interval_length

    def _generate_start_positon(self):
        rp = random.choice(list(self._graph.blocks.keys()))
        offset = random.randint(0, self._graph.node_size(rp)-1)
        return obg.Position(rp, offset)

    def generate_interval(self):
        start_pos = self._generate_start_positon()
        return self._generate_interval_from_start_pos(start_pos)

    def _generate_interval_from_start_pos(self, start_pos):
        cur_node = start_pos.region_path_id
        start_offset = start_pos.offset
        l = self._interval_length - (
            self._graph.node_size(cur_node)-start_offset)
        rps = [cur_node]
        while l > 0:
            cur_node = random.choice(self._graph.adj_list[cur_node])
            rps.append(cur_node)
            l -= self._graph.node_size(cur_node)

        end_offset = self._graph.node_size(cur_node) + l
        return obg.DirectedInterval(start_pos.offset, end_offset, rps)


class Path:
    def __init__(self, nodes, distances):
        self._nodes = np.array(nodes, dtype="int")
        self._distances = np.array(distances, dtype="int")
        self.length = distances[-1]

    def get_start_end(self, interval):
        region_paths = interval.region_paths
        if all(node in self._nodes for node in interval.region_paths):
            start = self._distances[list(self._nodes).index(region_paths[0])] + interval.start_position.offset
            end = self._distances[list(self._nodes).index(region_paths[-1])] + interval.end_position.offset
        elif all(-node in self._nodes for node in interval.region_paths):
            start = self._distances[list(self._nodes).index(-region_paths[-1])+1]-interval.end_position.offset
            end = self._distances[list(self._nodes).index(-region_paths[0])+1]-interval.start_position.offset
        else:
            return None

        return start, end

    def get_interval(self, start, end):
        is_reversed = start > end
        if is_reversed:
            start, end = end, start
        assert end <= self.length
        first_idx = np.where(self._distances > start)[0][0]-1
        last_idx = np.where(self._distances < end)[0][-1]
        rps = [self._nodes[i] for i in range(first_idx, last_idx+1)]
        if is_reversed:
            rps = [-rp for rp in rps[::-1]]
            start_offset = self._distances[last_idx+1]-end
            end_offset = self._distances[first_idx+1]-start
        else:
            start_offset = start-self._distances[first_idx]
            end_offset = end-self._distances[last_idx]
        return obg.DirectedInterval(start_offset, end_offset, rps)


class PathSimulator:
    def __init__(self, graph):
        self._graph = graph
        self._start_nodes = self._graph.get_first_blocks()

    def generate_path(self):
        start_node = random.choice(self._start_nodes)
        cur_node = start_node
        nodes = [start_node]
        sizes = [0, self._graph.node_size(start_node)]
        while self._graph.adj_list[cur_node]:
            cur_node = random.choice(
                self._graph.adj_list[cur_node])
            nodes.append(cur_node)
            sizes.append(self._graph.node_size(cur_node))
        dists = np.cumsum(sizes)
        print(len(nodes))
        return Path(nodes, dists)


class PeakReadSimulator:
    def __init__(self, genome_size, fragment_length, read_length):
        self._genome_size = genome_size
        self._fragment_length = fragment_length
        self._read_length = read_length
        self._summits = []

    def generate_summits(self, n_peaks=10):
        return (random.randint(self._fragment_length,
                               self._genome_size-self._fragment_length)
                for _ in range(n_peaks))

    def generate_read(self, summit):
        fragment_start = random.randint(
            summit - self._fragment_length, summit)
        if random.choice([True, False]):
            return fragment_start, fragment_start+self._fragment_length
        fragment_end = fragment_start + self._fragment_length
        return fragment_end, fragment_end-self._read_length

    def generate_reads_around_summit(self, summit, read_depth=100):
        self._summits.append(summit)
        return (self.generate_read(summit) for _ in range(read_depth))

    def generate_read_set(self, n_peaks=10, read_depth=100):
        return chain(*(self.generate_reads_around_summit(summit, read_depth)
                       for summit in self.generate_summits(n_peaks)))

    def get_summits(self):
        return self._summits


class EvaluateSimulations:
    data_folder = "../tests/"

    def __init__(self, graph_name="graph.obg"):
        self.graph = obg.GraphWithReversals.from_file(
            self.data_folder + graph_name)
        intervals = obg.IntervalCollection(
            list(self.generate_read_set(
                135, 36)))
        noise_generator = IntervalSimulator(self.graph, 36)
        noise = [noise_generator.generate_interval() for _ in range(2000)]
        intervals.intervals += noise
        print(len(intervals.intervals))
        for i in intervals:
            i.graph = self.graph
        info = ExperimentInfo(None, 135, 36)
        self.intervals = intervals.intervals
        self.map()
        callpeaks = CallPeaks(
            self.graph,
            vg_gam_file_to_interval_collection(None, "mapped_reads.gam", self.graph),  # obg.IntervalCollection((i for i in intervals)),
            vg_gam_file_to_interval_collection(None, "mapped_reads.gam", self.graph),  # obg.IntervalCollection((i for i in intervals)),
            info, has_control=False,
            out_file_base_name="simulated_",
            linear_map="../tests/haplo1kg50-mhc.lm")
        callpeaks.run()
        peaks = callpeaks.q_value_peak_caller.max_paths
        counter = 0
        for summit in self.summits:
            print(summit)
        for peak in peaks:
            if self.check_peak(peak):
                counter += 1
        print("%s/%s" % (counter, len(peaks)))

    def generate_read_set(self, fragment_length, read_length):
        path_generator = PathSimulator(self.graph)
        self.paths = [path_generator.generate_path() for _ in range(2)]
        simulators = [
            PeakReadSimulator(path.length, fragment_length, read_length)
            for path in self.paths]
        read_sets = [simulator.generate_read_set(80, 30)
                     for simulator in simulators]
        interval_sets = [[path.get_interval(*read) for read in read_set]
                         for path, read_set in zip(self.paths, read_sets)]
        self.summits = [simulator.get_summits() for simulator in simulators]
        intervals = interval_sets[0]+interval_sets[1]
        return intervals

    def check_peak(self, peak):
        ret = False
        letters = ["A", "B"]
        for path, summits, p in zip(self.paths, self.summits, letters):
            start_end = path.get_start_end(peak)
            if start_end is None:
                continue
            start, end = start_end
            covered_summits = [summit for summit in summits
                               if start < summit and end > summit]
            if covered_summits:
                print(p, ": ", peak, covered_summits[0])
                ret = True
        return ret

    def map(self, n=5):
        from pyvg.mapping import map
        self.write_sequence_and_intervals(n)
        graph = self.graph
        gam_file_name = "mapped_reads.gam"
        map("simulated_sequences.fq", "../tests/vgdata/haplo1kg50-mhc.xg",
            "../tests/vgdata/haplo1kg50-mhc.gcsa",
            gam_file_name)
        # reads_intervals = vg_gam_file_to_interval_collection(
        # None, gam_file_name, graph)

    def write_sequence_and_intervals(self, n=100):
        # obg.IntervalCollection(self.intervals).to_file(
        #     "simulated_intervals.py")
        logging.info("Getting sequences")
        self.retriever = SequenceRetriever.from_vg_graph(
            "../tests/haplo1kg50-mhc.vg")
        sequences = [self.retriever.get_interval_sequence(i)
                     for i in self.intervals]
        with open("simulated_sequences.fq", "w") as f:
            for i, seq in enumerate(sequences):
                f.write("@sim" + str(i) + "\n")
                f.write(seq + "\n")
                f.write("+\n")
                f.write("~"*36 + "\n")


class Comparer:
    def __init__(self):
        self.graph = obg.GraphWithReversals.from_file("../tests/graph.obg")

    def write_sequence_and_intervals(self, n=100):
        logging.info("Reading graph")
        logging.info("Simulating intervals")
        sim = IntervalSimulator(self.graph, 36)
        self.intervals = [sim.generate_interval() for _ in range(n)]
        obg.IntervalCollection(self.intervals).to_file(
            "simulated_intervals.py")
        logging.info("Getting sequences")
        self.retriever = SequenceRetriever.from_vg_graph(
            "../tests/haplo1kg50-mhc.vg")
        sequences = [self.retriever.get_interval_sequence(i)
                     for i in self.intervals]
        with open("simulated_sequences.fq", "w") as f:
            for i, seq in enumerate(sequences):
                f.write("@sim" + str(i) + "\n")
                f.write(seq + "\n")
                f.write("+\n")
                f.write("~"*36 + "\n")

    def compare_intervals(self, intervals_a, intervals_b):
        sequences_a = [self.retriever.get_interval_sequence(i)
                       for i in intervals_a]
        sequences_b = [self.retriever.get_interval_sequence(i)
                       for i in intervals_b]
        not_found = []
        print(len(sequences_a), len(sequences_b))
        for i, seq in enumerate(sequences_a):
            try:
                j = sequences_b.index(seq)
                sequences_b.pop(j)
                intervals_b.pop(j)
            except:
                print(seq)
                print(intervals_a[i])
                not_found.append(i)
        not_found = [sequences_a[i] for i in not_found]
        not_found.sort()
        print("---------------------------")
        print(len(not_found), len(sequences_b))
        for seq_a, seq_b in zip(not_found, sequences_b):
            print(seq_a)
            print(seq_b)

    def test_simulations(self, n=5):
        from pyvg.mapping import map
        from pyvg.util import vg_gam_file_to_interval_collection
        self.write_sequence_and_intervals(n)
        graph = obg.GraphWithReversals.from_file("../tests/graph.obg")
        gam_file_name = "mapped_reads.gam"
        map("simulated_sequences.fq", "../tests/vgdata/haplo1kg50-mhc.xg",
            "../tests/vgdata/haplo1kg50-mhc.gcsa",
            gam_file_name)
        reads_intervals = vg_gam_file_to_interval_collection(
             None, gam_file_name, graph)
        self.compare_intervals(list(self.intervals),
                               list(reads_intervals))

if __name__ == "__main__":
    random.seed(2000)
    EvaluateSimulations()
    exit()
    graph = obg.GraphWithReversals.from_file("../tests/graph.obg")
    intervals = obg.IntervalCollection(list(generate_read_set(
        graph,
        135, 36)))
    for i in intervals:
        i.graph = graph
    info = ExperimentInfo(None, 135, 36)
    callpeaks = CallPeaks(graph, intervals, intervals, info, has_control=False,
                          out_file_base_name="simulated_",
                          linear_map="../tests/haplo1kg50-mhc.lm")
    callpeaks.run()
