from itertools import chain
import numpy as np
import offsetbasedgraph as obg
from pyvg.sequences import SequenceRetriever
import random


class IntervalSimulator(object):
    def __init__(self, graph, interval_length):
        self._graph = graph
        self._interval_length = interval_length

    def _generate_start_positon(self):
        rp = random.choice(list(self._graph.blocks.keys()))
        offset = random.randint(0, self._graph.node_size(rp))
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
        return obg.Interval(start_pos.offset, end_offset, rps)


class Path:
    def __init__(self, nodes, distances):
        self._nodes = nodes
        self._distances = distances
        self.length = distances[-1]

    def get_interval(self, start, end):
        valid = np.logical_and(
            self._distances > start,
            self._distances < end)
        idxs = np.where(valid)[0]
        rps = self._nodes[idxs]
        start_offset = start-self._distances[idxs[0]-1]
        end_offset = end-self._distances[idxs[-1]]
        return obg.Interval(start_offset, end_offset, list(rps))


class PathSimulator:
    def __init__(self, graph):
        self._graph = graph
        self._start_nodes = self._graph.get_start_nodes()

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
        return Path(nodes, dists)


class PeakReadSimulator:
    def __init__(self, genome_size, fragment_length, read_length):
        self._genome_size = genome_size
        self._fragment_length = fragment_length
        self._read_length = read_length

    def generate_summits(self, n_peaks=10):
        return (random.randint(self._fragmen_length,
                               self._genome_size-self._fragment_length)
                for _ in range(n_peaks))

    def generate_read(self, summit):
        fragment_start = random.randint(
            summit - self._fragment_length, summit)
        if random.choice([True, False]):
            return fragment_start, fragment_start+self._length
        fragment_end = fragment_start + self._fragment_length
        return fragment_end, fragment_end-self._read_length

    def generate_reads_around_summit(self, summit, read_depth=100):
        return (self.generate_read(summit) for _ in range(read_depth))

    def generate_read_set(self, n_peaks=10, read_depth=100):
        return chain(*(self.generate_reads_around_summit(summit, read_depth)
                       for summit in self.generate_summits(n_peaks)))


def generate_read_set(graph, fragment_length, read_length):
    path_generator = PathSimulator(graph)
    paths = [path_generator.generate_path() for _ in range(2)]
    simulators = [PeakReadSimulator(path.length, fragment_length, read_length)
                  for path in paths]
    read_sets = [simulator.generate_read_set() for simulator in simulators()]
    interval_sets = [(path.get_interval(*read) for read in read_set)
                     for path, read_set in zip(paths, read_sets)]
    return chain(*interval_sets)


if __name__ == "__main__":
    graph = obg.GraphWithReversals.from_file("../tests/graph.obg")
    sim = IntervalSimulator(graph, 36)
    intervals = [sim.generate_interval() for _ in range(100)]
    obg.IntervalCollection(intervals).to_file("simulated_intervals.py")
    retriever = SequenceRetriever.from_vg_graph("../tests/haplo1kg50-mhc.vg")
    sequences = [retriever.get_interval_sequence(i) for i in intervals]
    with open("simulated_sequences.fq", "w") as f:
        for i, seq in enumerate(sequences):
            f.write("@sim" + str(i) + "\n")
            f.write(seq + "\n")
            f.write("+\n")
            f.write("~"*36 + "\n")
