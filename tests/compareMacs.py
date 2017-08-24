from offsetbasedgraph import Graph, Block, Position, Interval
from offsetbasedgraph.interval import IntervalCollection
from graph_peak_caller.callpeaks import CallPeaks, ExperimentInfo
import subprocess
import random
import re
import numpy as np


class SimpleInterval(object):

    def __init__(self, start, end, direction):
        self.node_id = 0
        self.start = start
        self.end = end
        self.direction = direction
        if direction is not None:
            self.dir_symbol = "+" if direction > 0 else "-"

    def to_file_line(self):
        return "\t".join(
            str(v) for v in
            (self.node_id, self.start, self.end, ".", "0", self.dir_symbol))+"\n"

    @classmethod
    def from_file_line(cls, line):
        _, start, end, _, _, dir_symbol = line.split("\t")
        direction = +1 if dir_symbol == "+" else -1
        return cls(int(start), int(end), direction)

    def __eq__(self, other):
        attrs = ["start", "end", "direction"]
        return (getattr(self, attr) == getattr(other, attr)
                for attr in attrs)


class ValuedInterval(SimpleInterval):
    def __init__(self, start, end, value):
        SimpleInterval.__init__(self, start, end, None)
        self.value = value

    @classmethod
    def from_file_line(cls, line):
        node_id, start, end, value = line.split("\t")
        obj = cls(int(start), int(end), value)
        obj.node_id = int(node_id)
        return obj


class MACSTests(object):
    def __init__(self, node_size, n_nodes, n_intervals, read_length=15):
        self.node_size = node_size
        self.n_nodes = n_nodes
        self.n_intervals = n_intervals
        self.read_length = read_length
        self.genome_size = node_size*n_nodes
        self.setup()

    def setup(self):
        self.create_linear_graph()
        self.create_intervals()
        self.write_intervals()

    def linear_to_graph_interval(self, lin_interval):
        start = lin_interval.start
        end = lin_interval.end
        start_rp = start//self.node_size
        end_rp = (end-1)//self.node_size
        start_pos = Position(
            start_rp,
            start % self.node_size)
        end_pos = Position(
            end_rp,
            ((end-1) % self.node_size) + 1)
        region_paths = list(range(start_rp, end_rp+1))
        return Interval(start_pos, end_pos, region_paths,
                        direction=lin_interval.direction, graph=self.graph)

    def _convert_valued_interval(self, interval):
        interval.start += self.node_size*interval.node_id
        interval.end += self.node_size*interval.node_id

    def graph_to_linear_pos(self, pos):
        return pos.region_path_id*self.node_size + pos.offset

    def graph_to_linear_interval(self, graph_interval):
        start = self.graph_to_linear_pos(graph_interval.start_position)
        end = self.graph_to_linear_pos(graph_interval.end_position)
        return SimpleInterval(start, end, graph_interval.direction)

    def assertEqualIntervals(self, linear_intervals, graph_intervals):
        graph_intervals = [self.graph_to_linear_interval(g_interval)
                           for g_interval in graph_intervals]
        assert len(graph_intervals) == len(linear_intervals)
        for interval in graph_intervals:
            assert interval in linear_intervals

    def assertEqualIntervalFiles(self, graph_file, linear_file):
        graph_intervals = IntervalCollection.create_generator_from_file(
            graph_file)
        linear_intervals = (SimpleInterval.from_file_line(line) for
                            line in open(linear_file).readlines())
        self.assertEqualIntervals(list(linear_intervals), graph_intervals)

    def test_filter_dup(self):
        caller = CallPeaks("lin_graph", "graph_intervals")
        command = "macs2 filterdup -i %s --keep-dup=1 -o %s" % (
            "lin_intervals.bed", "lin_intervals_dup.bed")
        print(command)
        command = command.split()
        subprocess.check_output(command)
        self.dup_file_name = caller.filter_duplicates("graph_intervals")
        self.assertEqualIntervalFiles(
            self.dup_file_name,
            "lin_intervals_dup.bed")

    def _create_pileup(self, pileup_file, convert=False):
        pileup = np.zeros(self.genome_size)
        valued_intervals = (ValuedInterval.from_file_line(line) for line in
                            open(pileup_file).readlines())
        for interval in valued_intervals:
            if convert:
                self._convert_valued_interval(interval)
            pileup[interval.start:interval.end] = interval.value
        return pileup

    def assertPileupFilesEqual(self, graph_file, linear_file):
        linear_pileup = self._create_pileup(linear_file)
        graph_pileup = self._create_pileup(graph_file, convert=True)
        print(linear_pileup)
        print(graph_pileup)
        assert all(linear_pileup == graph_pileup)

    def test_sample_pileup(self):
        info = ExperimentInfo(self.genome_size, self.n_intervals,
                              self.n_intervals, 100, self.read_length)
        caller = CallPeaks("lin_graph", self.dup_file_name,
                           experiment_info=info)
        caller.create_graph()
        caller.preprocess()
        caller.create_sample_pileup()
        command = "macs2 pileup -i %s -o %s --extsize %s" % (
            "lin_intervals_dup.bed", "lin_sample_pileup.bdg", 100)
        print(command)
        subprocess.check_output(command.split())
        print(caller._sample_track)
        self.assertPileupFilesEqual(
            caller._sample_track,
            "lin_sample_pileup.bdg")

    def write_intervals(self):
        f = open("lin_intervals.bed", "w")
        f.writelines(interval.to_file_line() for
                     interval in self.linear_intervals)
        f.close()
        graph_intervals = IntervalCollection(self.graph_intervals)
        graph_intervals.to_file("graph_intervals")

    def create_linear_graph(self):
        nodes = {i: Block(self.node_size) for i in range(self.n_nodes)}
        adj_list = {i: [i+1] for i in range(self.n_nodes-1)}
        self.graph = Graph(nodes, adj_list)
        self.graph.to_file("lin_graph")

    def _get_graph_interval(self, start, end, direction):
        start_rp = start//self.node_size
        end_rp = (end+1)//self.node_size
        start_pos = Position(
            start_rp,
            start % self.node_size)
        end_pos = Position(
            end_rp,
            (end % self.node_size) + 1)
        region_paths = list(range(start_rp, end_rp+1))
        return Interval(start_pos, end_pos, region_paths, direction=direction)

    def create_intervals(self):
        self.linear_intervals = []
        self.graph_intervals = []
        for _ in range(self.n_intervals):
            direction = random.choice((-1, 1))
            start = random.randint(0, self.genome_size-self.read_length)
            end = start+self.read_length
            interval = SimpleInterval(start, end, direction)
            self.linear_intervals.append(interval)
            self.graph_intervals.append(self.linear_to_graph_interval(interval))
            if start % 10 == 0:
                self.linear_intervals.append(interval)
                self.graph_intervals.append(self.linear_to_graph_interval(interval))

    def test_shift_estimation(self):
        self.setup()
        caller = CallPeaks("lin_graph", "graph_intervals")
        caller.create_graph()
        caller.find_info()
        read_length_graph = caller.read_length
        fragment_length_graph = caller.fragment_length

        # Macs
        command = ["macs", "predictd", "-i", "lin_intervals.bed", "-g",
                   self.genome_size, "m", "5", "50"]
        output = subprocess.check_output(command, stderr=subprocess.STDOUT)
        output = output.decode("utf-8")
        print(output)
        tag_size = re.search("tag size = ([0-9]+)", output).groups()[0]
        tag_size = int(tag_size)
        fragment_length = re.search("fragment length is ([0-9]+) bp", output).groups()[0]
        fragment_length = int(fragment_length)

        assert read_length_graph == tag_size
        assert fragment_length_graph == fragment_length

if __name__ == "__main__":
    test = MACSTests(100, 100, 100)
    test.test_filter_dup()
    test.test_sample_pileup()
    # test.test_shift_estimation()
