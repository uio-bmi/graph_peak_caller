import cProfile
import pstats
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

    def __str__(self):
        return "(%s-%s:%s)" % (self.node_id, self.start, self.end)

    __repr__ = __str__

    def to_file_line(self):
        return "\t".join(
            str(v) for v in
            (self.node_id, self.start, self.end, ".", "0", self.dir_symbol))+"\n"

    @classmethod
    def from_file_line(cls, line):
        if line.startswith("track"):
            return None
        node_id, start, end, _, _, dir_symbol = line.split("\t")[:6]
        direction = +1 if dir_symbol == "+" else -1
        if dir_symbol == ".":
            direction = None

        obj = cls(int(start), int(end), direction)
        obj.node_id = int(node_id)
        return obj

    def __eq__(self, other):
        attrs = ["start", "end", "direction"]
        return (getattr(self, attr) == getattr(other, attr)
                for attr in attrs)
    def __str__(self):
        return "%s %d %d %s" % (self.node_id, self.start, self.end, self.dir_symbol)

    def __repr__(self):
        return self.__str__()


class ValuedInterval(SimpleInterval):
    def __init__(self, start, end, value):
        SimpleInterval.__init__(self, start, end, None)
        self.value = value

    @classmethod
    def from_file_line(cls, line):
        if line.startswith("track"):
            return None
        node_id, start, end, value = line.split("\t")
        obj = cls(int(start), int(end), value)
        obj.node_id = int(node_id)
        return obj


class MACSTests(object):
    def __init__(self, node_size, n_nodes, n_intervals, read_length=15, fragment_length=50):
        self.node_size = node_size
        self.n_nodes = n_nodes
        self.n_intervals = n_intervals
        self.read_length = read_length
        self.genome_size = node_size*n_nodes
        self.fragment_length = fragment_length
        self.setup()

    def setup(self):
        self.create_linear_graph()
        self.create_intervals()
        self.write_intervals()
        self.info = ExperimentInfo(self.genome_size, self.n_intervals,
                                   self.n_intervals, self.fragment_length,
                                   self.read_length)
        self.caller = CallPeaks("lin_graph", "graph_intervals",
                                experiment_info=self.info,
                                verbose=True)
        self.caller.create_graph()

    # Tests
    def test_filter_dup(self):
        command = "macs2 filterdup -i %s --keep-dup=1 -o %s" % (
            "lin_intervals.bed", "lin_intervals_dup.bed")
        command = command.split()
        subprocess.check_output(command)
        self.dup_file_name = self.caller.filter_duplicates("graph_intervals", write_to_file="graph_intervals_filtered")
        self.assertEqualIntervalFiles(
            self.dup_file_name,
            "lin_intervals_dup.bed")

    def test_sample_pileup(self):
        # self.caller.create_graph()
        self.caller.create_sample_pileup(True)
        self._create_sample_pileup()
        self.assertPileupFilesEqual(
            self.caller._sample_track,
            "lin_sample_pileup.bdg")

    def test_control_pileup(self):
        self.caller.create_control(True)
        self._create_control()
        assert isinstance(self.caller._control_track, str)
        self.assertPileupFilesEqual(self.caller._control_track,
                                    "lin_control_pileup.bdg", min_value=self.background)

    def test_call_peaks(self):
        self.caller.get_p_values()
        self._get_scores()
        # self.assertPileupFilesEqual(self.caller._p_value_track,
        # "lin_scores.bdg")

        self.caller.call_peaks()
        self._call_peaks()
        self.assertEqualBedFiles("final_peaks", "lin_peaks.bed")

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
        assert len(graph_intervals) == len(linear_intervals), \
                "%d != %d" % (len(graph_intervals), len(linear_intervals))
        for interval in graph_intervals:
            assert interval in linear_intervals

    def assertEqualIntervalFiles(self, graph_file, linear_file):
        graph_intervals = IntervalCollection.create_generator_from_file(
            graph_file)
        linear_intervals = (SimpleInterval.from_file_line(line) for
                            line in open(linear_file).readlines())
        self.assertEqualIntervals(list(linear_intervals), graph_intervals)

    def _create_binary_track(self, intervals):
        pileup = np.zeros(self.genome_size, dtype="bool")
        for interval in intervals:
            if interval is None:
                continue
            pileup[interval.start:interval.end] = True
        return pileup

    def assertEqualBedFiles(self, graph_file, linear_file):
        graph_intervals = [SimpleInterval.from_file_line(line)
                           for line in open(graph_file).readlines()]
        linear_intervals = [SimpleInterval.from_file_line(line) for
                            line in open(linear_file).readlines()]

        for graph_interval in graph_intervals:
            self._convert_valued_interval(graph_interval)
        pileup1 = self._create_binary_track(linear_intervals)
        pileup2 = self._create_binary_track(graph_intervals)
        indices = np.where(pileup1 != pileup2)
        #for i in indices[0]:
        #    print(i)

        print(indices)
        print("Pileup1")
        print(pileup1)
        print("Pileup2")
        print(pileup2)
        assert np.allclose(pileup1, pileup2)

    def _create_pileup(self, pileup_file, convert=False, limit=False, min_value=None):
        pileup = np.zeros(self.genome_size)
        valued_intervals = (ValuedInterval.from_file_line(line) for line in
                            open(pileup_file).readlines())
        for interval in valued_intervals:
            if interval is None:
                continue
            if convert:
                self._convert_valued_interval(interval)
            pileup[interval.start:interval.end] = interval.value
        if min_value is not None:
            pileup = np.maximum(pileup, min_value)
        return pileup

    def assertPileupFilesEqual(self, graph_file, linear_file, min_value=None):
        assert isinstance(graph_file, str)
        assert isinstance(linear_file, str)

        linear_pileup = self._create_pileup(linear_file, min_value=min_value)
        graph_pileup = self._create_pileup(graph_file, convert=True)
        # assert not all(graph_pileup == graph_pileup[0])
        assert sum(graph_pileup) > 0

        if not np.allclose(linear_pileup, graph_pileup):
            print(linear_pileup)
            print(graph_pileup)
            print(np.where(linear_pileup != graph_pileup))
        assert np.allclose(linear_pileup, graph_pileup), \
            "Pileup in %s != pileup in %s" % (linear_file, graph_file)

    def _create_sample_pileup(self):
        command = "macs2 pileup -i %s -o %s --extsize %s" % (
            "lin_intervals.bed", "lin_sample_pileup.bdg", self.fragment_length)
        print(command)
        subprocess.check_output(command.split())

    def _get_scores(self):
        command = "macs2 bdgcmp -t lin_sample_pileup.bdg -c lin_control_pileup.bdg -m ppois -o lin_scores.bdg"
        print(command)
        subprocess.check_output(command.split())

    def _call_peaks(self):
        threshold = -np.log10(0.05)
        command = "macs2 bdgpeakcall -i lin_scores.bdg -c %s -l %s -g %s -o lin_peaks.bed" % (
            threshold, self.fragment_length, self.read_length)
        print(command)
        subprocess.check_output(command.split())


    def _create_control(self):
        for ext in [2500, 5000]:
            command = "macs2 pileup -i %s -o %s -B --extsize %s" % (
                "lin_intervals.bed", "lin_control_pileup%s.bdg" % ext, ext)
            subprocess.check_output(command.split())
            command = "macs2 bdgopt -i lin_control_pileup%s.bdg -m multiply -p %s -o lin_control_pileup%s.bdg" % (
                ext, self.fragment_length/(ext*2), ext)
            subprocess.check_output(command.split())
        command = "macs2 bdgcmp -m max -t lin_control_pileup2500.bdg -c lin_control_pileup5000.bdg -o lin_control_pileup.bdg"

        subprocess.check_output(command.split())

        self.background = self.n_intervals * self.fragment_length / self.genome_size
        command = "macs2 bdgopt -i lin_control_pileup.bdg -m max -p %s -o lin_control_pileup.bdg" % self.background
        subprocess.check_output(command.split())

    def write_intervals(self):
        f = open("lin_intervals.bed", "w")
        f.writelines(interval.to_file_line() for
                     interval in self.linear_intervals)
        f.close()
        print("Wrote to lin_intervals.bed")
        graph_intervals = IntervalCollection(self.graph_intervals)
        graph_intervals.to_file("graph_intervals")
        print("Wrote to graph_intervals")

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
        self.n_duplicates = 0
        i = 0
        random.seed(1)
        for _ in range(self.n_intervals):
            direction = random.choice((-1, 1))
            if direction == -1:
                start = random.randint(self.fragment_length-self.read_length,
                                       self.genome_size-self.read_length)
            else:
                start = random.randint(0, self.genome_size-self.read_length-self.fragment_length)
            end = start+self.read_length
            interval = SimpleInterval(start, end, direction)
            self.linear_intervals.append(interval)
            self.graph_intervals.append(self.linear_to_graph_interval(interval))

            i += 1

            # Add duplicate
            if start % 10 == 0:
                self.linear_intervals.append(interval)
                self.graph_intervals.append(self.linear_to_graph_interval(interval))
                self.n_duplicates += 1
                i += 1

            # Add pair
            if (start % 2 == 0) or True:
                if direction == 1:
                    start = start + self.fragment_length - self.read_length
                else:
                    start = end - self.fragment_length

                end = start + self.read_length
                direction = direction * -1
                paired_interval = SimpleInterval(start, end, direction)
                self.linear_intervals.append(paired_interval)
                self.graph_intervals.append(self.linear_to_graph_interval(paired_interval))
                #print("Creating paired interval")
                #print(paired_interval)
                #print("to")
                #print(interval)
                i += 1

        self.n_intervals = i
        print("Created totalt %d intervals " % i)

        assert len(self.linear_intervals) == len(self.graph_intervals)

    def test_shift_estimation(self):
        self.setup()
        caller = CallPeaks("lin_graph", "graph_intervals")
        caller.create_graph()
        info = ExperimentInfo.find_info(caller.ob_graph, caller.sample_file_name, caller.control_file_name)
        read_length_graph = info.read_length
        fragment_length_graph = info.fragment_length

        # Macs
        command = ["macs2", "predictd", "-i", "lin_intervals.bed", "-g", str(self.genome_size), "-m", "5", "50"]
        string_commmand = ' '.join(command)
        print(string_commmand)
        output = subprocess.check_output(command, stderr=subprocess.STDOUT)
        output = output.decode("utf-8")
        print(output)
        tag_size = re.search("tag size = ([0-9]+)", output).groups()[0]
        tag_size = int(tag_size)
        fragment_length = re.search("fragment length is ([0-9]+) bp", output).groups()[0]
        fragment_length = int(fragment_length)

        assert read_length_graph == tag_size, \
            "Read length from graph % d != %d (macs reads length)" % (read_length_graph, tag_size)
        assert fragment_length_graph == fragment_length

    def profile(self):
        self.caller.run()

    def _run_whole_macs(self):
        command = "macs2 callpeak -t lin_intervals.bed -f BED -g " + str(self.genome_size) + " -n macstest -B -p 0.05 --slocal 5000 --llocal 10000"
        print(command)
        command = command.split()
        output = subprocess.check_output(command, stderr=subprocess.STDOUT)
        output = output.decode("utf-8")
        print(output)

    def test_whole_pipeline(self):
        #self.caller.run("final_peaks")
        self.caller.create_graph()

        self.caller.info = ExperimentInfo.find_info(
            self.caller.ob_graph, self.caller.sample_file_name, self.caller.control_file_name)

        self.caller.preprocess()

        self.caller.create_sample_pileup()

        self.assertPileupFilesEqual("sample_track.bdg", "macstest_treat_pileup.bdg")
        self.caller.create_control(True)
        #self.caller.scale_tracks()
        #self.caller.get_p_values()
        #self.caller.call_peaks("final_peaks")


        #self._run_whole_macs()
        #self.assertPileupFilesEqual("control_track.bdg", "macstest_control_lambda.bdg")
        #self.assertEqualBedFiles("final_peaks", "test_peaks.narrowPeak")

def small_test():
    return MACSTests(1000, 1000, 100000, read_length=15, fragment_length=20)


def big_test():
    return MACSTests(5000, 10000, 100000, read_length=51, fragment_length=120)


if __name__ == "__main__":
    test = big_test()
    #cProfile.run("test.profile()", "profiling")
    #p = pstats.Stats("profiling")
    #p.sort_stats("tottime").print_stats()
    #exit()

    #test.test_filter_dup()
    #test.test_shift_estimation()
    #test.test_sample_pileup()
    #test.test_control_pileup()
    #test.test_call_peaks()
    test.test_whole_pipeline()