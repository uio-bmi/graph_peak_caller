import cProfile
import pstats
import subprocess
import random
import re
import numpy as np
import logging

from offsetbasedgraph import Graph, Block, Position,\
    DirectedInterval
from offsetbasedgraph.interval import IntervalCollection
from graph_peak_caller.callpeaks import CallPeaks, ExperimentInfo
from graph_peak_caller.pileup import Pileup
from graph_peak_caller.sparsepileup import SparsePileup

logging.basicConfig(level=logging.ERROR)


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
            (self.node_id, self.start, self.end, ".",
             "0", self.dir_symbol)) + "\n"

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
        return "%s %d %d %s" % (
            self.node_id, self.start, self.end, self.dir_symbol)

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
    def __init__(self, node_size, n_nodes, n_intervals,
                 read_length=15, fragment_length=50, with_control=False):
        self.node_size = node_size
        self.n_nodes = n_nodes
        self.with_control = with_control
        self.n_intervals = n_intervals
        self.read_length = read_length
        self.genome_size = node_size*n_nodes
        self.fragment_length = fragment_length
        self.setup()

    def setup(self):
        print("######## SETUP ########")
        self.create_linear_graph()
        self.create_intervals()
        self.write_intervals()
        self.n_intervals_control = self.n_intervals
        self.info = ExperimentInfo(self.genome_size,
                                   self.fragment_length - 1 - (self.fragment_length%2),
                                   self.read_length)
        self.info.n_control_reads = self.n_intervals
        self.info.n_sample_reads = self.n_intervals

        control_file_name = "graph_intervals"
        if self.with_control:
            control_file_name = "graph_intervals_control"

        self.caller = CallPeaks("lin_graph", "graph_intervals", control_file_name,
                                has_control=self.with_control,
                                experiment_info=self.info,
                                verbose=True)

        self.caller.create_graph()

    # Tests
    def test_filter_dup(self):
        command = "macs2 filterdup -i %s --keep-dup=1 -o %s" % (
            "lin_intervals.bed", "lin_intervals_dup.bed")
        command = command.split()
        subprocess.check_output(command)
        self.dup_file_name = self.caller.filter_duplicates(
            "graph_intervals",
            write_to_file="graph_intervals_filtered")
        self.assertEqualIntervalFiles(
            self.dup_file_name,
            "lin_intervals_dup.bed")

    def test_sample_pileup(self):
        # self.caller.create_graph()
        self.caller.sample_intervals = self.graph_intervals
        self.caller.create_sample_pileup(True)
        self._create_sample_pileup()
        self.assertPileupFilesEqual(
            self.caller._sample_track,
            "lin_sample_pileup.bdg")

    def test_control_pileup(self):
        self.caller.control_intervals = self.graph_intervals
        self.caller.create_control(True)
        self._create_control()
        assert isinstance(self.caller._control_track, str)
        self.assertPileupFilesEqual(self.caller._control_track,
                                    "lin_control_pileup.bdg", min_value=self.background)

    def test_call_peaks(self):
        # self.assertPileupFilesEqual("control_track.bdg", "macstest_control_lambda.bdg")
        # self.assertPileupFilesEqual("sample_track.bdg", "macstest_treat_pileup.bdg")
        self.caller._control_pileup = SparsePileup.from_bed_graph(self.graph, "control_track.bdg")
        self.caller._sample_pileup = SparsePileup.from_bed_graph(self.graph, "sample_track.bdg")
        self.caller.get_score()
        # self.caller.p_values.to_bed_graph(self.caller._p_value_track)
        self._get_scores("qpois")
        #self.assertPileupFilesEqual(self.caller._p_value_track,
        # "lin_scores.bdg")

        # self.caller.get_q_values()
        self.caller.q_values.to_bed_graph(self.caller._q_value_track)
        self._get_scores()
        print(self.caller._q_value_track)
        self.assertPileupFilesEqual(self.caller._q_value_track,
                                    "lin_scores.bdg")
        self._call_peaks()
        self.caller.call_peaks()

        self.assertEqualBedFiles("final_peaks", "lin_peaks.bed")

    def neg_linear_to_graph_interval(self, lin_interval):
        start_offset = (-lin_interval.end) % self.node_size
        end_offset = (-lin_interval.start+1) % self.node_size - 1
        start_rp = (lin_interval.end-1) // self.node_size + 1
        end_rp = (lin_interval.start) // self.node_size + 1
        rps = list(range(start_rp*-1, end_rp*-1+1))
        return DirectedInterval(start_offset, end_offset, rps,
                                graph=self.graph)

    def linear_to_graph_interval(self, lin_interval):
        if lin_interval.direction == -1:
            return self.neg_linear_to_graph_interval(lin_interval)
        start = lin_interval.start
        end = lin_interval.end
        start_rp = start//self.node_size+1
        end_rp = (end-1)//self.node_size+1
        start_pos = Position(
            start_rp,
            start % self.node_size)
        end_pos = Position(
            end_rp,
            ((end-1) % self.node_size) + 1)
        region_paths = list(range(start_rp, end_rp+1))
        interval = DirectedInterval(
            start_pos, end_pos, region_paths,
            direction=lin_interval.direction, graph=self.graph)
        return interval

    def _convert_valued_interval(self, interval):
        true_id = abs(interval.node_id)-1
        interval.start += self.node_size*true_id
        interval.end += self.node_size*true_id

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
        graph_intervals = IntervalCollection.from_file(
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

        print("Pileup1")
        print(pileup1[indices])
        print("Pileup2")
        print(pileup2[indices])
        assert np.allclose(pileup1, pileup2)

    def _create_pileup(self, pileup_file, convert=False, limit=False,
                       min_value=None):
        pileup = np.zeros(self.genome_size)
        valued_intervals = (ValuedInterval.from_file_line(line) for line in
                            open(pileup_file).readlines())
        for interval in valued_intervals:
            if interval is None:
                continue
            if convert:
                self._convert_valued_interval(interval)
            try:
                pileup[interval.start:interval.end] = interval.value
            except:
                print(interval.start, interval.end, interval.value)
                raise
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
        rtol = 0.001

        if not np.allclose(linear_pileup, graph_pileup, rtol=rtol):
            different = np.abs(linear_pileup - graph_pileup) > rtol
            print(linear_pileup)
            print(graph_pileup)
            print(different)
            print(np.where(different))
            print("Number of indices different")
            print(len(np.where(different)[0]))
            if not len(np.where(different)[0]):
                return
            print("Differences:")

            #for different_index in np.where(different)[0]:
            #    print("%.3f %.3f" % (linear_pileup[different_index], graph_pileup[different_index]))


            #print(linear_pileup[np.where(different)[0]])
            #print(graph_pileup[np.where(different)[0]])

        assert np.allclose(linear_pileup[:-250], graph_pileup[:-250], rtol=rtol), \
            "Pileup in %s != pileup in %s" % (linear_file, graph_file)

    def _create_sample_pileup(self):
        command = "macs2 pileup -i %s -o %s --extsize %s -f BED" % (
            "lin_intervals.bed", "lin_sample_pileup.bdg", self.fragment_length - 1)
        print(command)
        subprocess.check_output(command.split())

    def _get_scores(self, t="qpois"):
        command = "macs2 bdgcmp -t lin_sample_pileup.bdg -c lin_control_pileup.bdg -m %s -o lin_scores.bdg" % t
        # command = "macs2 bdgcmp -t macstest_treat_pileup.bdg -c macstest_control_lambda.bdg  -m %s -o lin_scores.bdg" % t
        print(command)
        subprocess.check_output(command.split())

    def _call_peaks(self):
        threshold = -np.log10(0.05)
        command = "macs2 bdgpeakcall -i lin_scores.bdg -c %s -l %s -g %s -o lin_peaks.bed" % (
            threshold, self.info.fragment_length, self.read_length)
        print(command)
        subprocess.check_output(command.split())

    def _create_control(self):
        for ext in [2500]:
            command = "macs2 pileup -i %s -o %s -B --extsize %s" % (
                "lin_intervals.bed", "lin_control_pileup%s.bdg -f BED" % ext, ext)
            subprocess.check_output(command.split())
            command = "macs2 bdgopt -i lin_control_pileup%s.bdg -m multiply -p %s -o lin_control_pileup%s.bdg" % (
                ext, (self.fragment_length-1)/(ext*2), ext)
            subprocess.check_output(command.split())
        # command = "macs2 bdgcmp -m max -t lin_control_pileup2500.bdg -c lin_control_pileup5000.bdg -o lin_control_pileup.bdg"

        # subprocess.check_output(command.split())

        self.background = self.n_intervals * self.info.fragment_length / self.genome_size
        logging.info(self.background)
        command = "macs2 bdgopt -i lin_control_pileup2500.bdg -m max -p %s -o lin_control_pileup.bdg" % self.background
        print(command)
        subprocess.check_output(command.split())

    def write_intervals(self):
        f = open("lin_intervals.bed", "w")
        f.writelines(interval.to_file_line() for
                     interval in self.linear_intervals)
        f.close()
        print("Wrote to lin_intervals.bed")
        graph_intervals = IntervalCollection(self.graph_intervals)
        graph_intervals.to_file("graph_intervals", True)

        if self.with_control:
            f = open("lin_intervals_control.bed", "w")
            f.writelines(interval.to_file_line() for
                         interval in self.linear_intervals_control)
            f.close()
            graph_intervals = IntervalCollection(self.graph_intervals_control)
            graph_intervals.to_file("graph_intervals_control")

        print("Wrote to graph_intervals")

    def create_linear_graph(self):
        nodes = {i+1: Block(self.node_size) for i in range(self.n_nodes)}
        adj_list = {i: [i+1] for i in range(1, self.n_nodes)}
        self.graph = Graph(nodes, adj_list)
        self.graph.to_file("lin_graph")

    def _get_graph_interval(self, tmp_start, tmp_end, direction):
        start = tmp_start
        end = tmp_end
        if direction == -1:
            start = -tmp_end
            end = -tmp_start
        start_rp = start//self.node_size
        end_rp = (end+1)//self.node_size
        region_paths = list(range(start_rp, end_rp))
        start_pos = Position(
            start_rp,
            start % self.node_size)
        end_pos = Position(
            end_rp,
            (end % self.node_size) + 1)
        return DirectedInterval(start_pos, end_pos, region_paths, direction=direction)

    def create_pairs_around_point(self, point, n=1):
        intervals = []
        for _ in range(n):
            offset = random.randint(-n, n)
            point = point+offset
            pos_start = point-self.fragment_length//2
            pos_end = pos_start+self.read_length
            if pos_start > 0 and pos_end < self.genome_size:
                intervals.append(SimpleInterval(pos_start, pos_end, 1))
                assert pos_start >= 0 and pos_end >= 0
            neg_end = point+self.fragment_length//2
            neg_start = neg_end-self.read_length
            if neg_end < self.genome_size and neg_start >= 0:
                intervals.append(SimpleInterval(neg_start, neg_end, -1))
                assert neg_start >= 0 and neg_end >= 0

        return intervals

    def create_random_linear_reads(self, n_reads, include_pairs=False):
        reads = []
        for i in range(n_reads//100+1):
            point = random.randint(0, self.genome_size)
            reads.extend(self.create_pairs_around_point(point))
            continue
            direction = random.choice((-1, 1))
            if direction == -1:
                start = random.randint(self.fragment_length-self.read_length,
                                       self.genome_size-self.read_length)
            else:
                start = random.randint(0, self.genome_size-self.read_length-self.fragment_length)
            end = start+self.read_length
            interval = SimpleInterval(start, end, direction)
            reads.append(interval)

            # Add duplicate
            if start % 50 == 0:
                reads.append(interval)

            # Add pair
            if include_pairs:
                if direction == 1:
                    start = start + self.fragment_length - self.read_length
                else:
                    start = end - self.fragment_length

                end = start + self.read_length
                direction = direction * -1
                paired_interval = SimpleInterval(start, end, direction)
                reads.append(paired_interval)

        return reads

    def create_intervals(self):
        self.linear_intervals = self.create_random_linear_reads(self.n_intervals, include_pairs=True)
        dummy_end = SimpleInterval(self.genome_size - self.read_length, self.genome_size, -1)
        self.linear_intervals.append(dummy_end)
        self.graph_intervals = [self.linear_to_graph_interval(i) for i in self.linear_intervals]
        logging.debug(len(self.graph_intervals))
        # assert all(i.start >= 0 and i.end >= 0 for i in self.linear_intervals)
        # t = all(all(abs(rp) in self.graph.blocks for rp in read.region_paths) for read in self.graph_intervals)
        # print([i for i in self.graph_intervals if not all(abs(rp) in self.graph.blocks for rp in i.region_paths)])
        # assert t

        self.n_intervals = len(self.linear_intervals)
        self.linear_intervals = sorted(self.linear_intervals, key = lambda x: (x.node_id, x.start))
        self.graph_intervals = sorted(self.graph_intervals, key = lambda x: (x.region_paths[0], x.start_position.offset))
        print("Created %d intervals " % self.n_intervals)

        if self.with_control:
            self.linear_intervals_control = self.create_random_linear_reads(self.n_intervals // 2, include_pairs=False)

            self.linear_intervals_control.append(dummy_end)

            self.graph_intervals_control = [self.linear_to_graph_interval(i) for i in self.linear_intervals_control]
            self.n_intervals_control = len(self.linear_intervals_control)
            print("Created %d control intervals " % self.n_intervals_control)
        else:
            self.n_intervals_control = self.n_intervals

    def test_shift_estimation(self):
        self.setup()
        caller = CallPeaks("lin_graph", "graph_intervals_filtered", "graph_intervals_filtered", has_control=False)
        caller.create_graph()
        info = ExperimentInfo.find_info(
            caller.ob_graph, caller.sample_file_name, caller.control_file_name)
        read_length_graph = info.read_length
        fragment_length_graph = info.fragment_length

        # Macs
        command = ["macs2", "predictd", "-i", "lin_intervals_dup.bed", "-g", str(self.genome_size), "-m", "5", "50"]
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

        command = "macs2 callpeak -t lin_intervals.bed -f BED -g " + str(self.genome_size) + " -n macstest -B -q 0.05 --llocal 5000"
        if self.with_control:
            command += " --slocal 2500 -c lin_intervals_control.bed"

        print(command)
        command = command.split()
        output = subprocess.check_output(command, stderr=subprocess.STDOUT)
        output = output.decode("utf-8")
        print(output)

    def test_whole_pipeline(self):
        self._run_whole_macs()
        self.caller.create_graph()
        self.caller.preprocess()

        self.caller.create_sample_pileup(True)
        print("#", self.caller.info.n_control_reads)
        self.caller.create_control(True)
        print("#", self.caller.info.n_control_reads)
        self.caller.scale_tracks(update_saved_files=True)
        self.assertPileupFilesEqual("sample_track.bdg", "macstest_treat_pileup.bdg")
        self.assertPileupFilesEqual("control_track.bdg", "macstest_control_lambda.bdg")
        print("################### GETTING SCORE")
        self.caller.get_score()
        print("################### CALLING PEAKS")
        self.caller.call_peaks("final_peaks")

        self.assertEqualBedFiles("final_peaks", "macstest_peaks.narrowPeak")

    def test_final_tracks(self):
        self._run_whole_macs()
        self.caller.run()
        self.assertEqualBedFiles("final_peaks", "macstest_peaks.narrowPeak")


def small_test(with_control=False):
    return MACSTests(1000, 10, 1000, read_length=10,
                     fragment_length=30, with_control=with_control)


def big_test(with_control=False):
    return MACSTests(100, 1000, 100000, read_length=51,
                     fragment_length=120, with_control=with_control)


if __name__ == "__main__":
    random.seed(102)
    test = big_test(False)
    test.test_sample_pileup()
    test.test_control_pileup()
    test.test_call_peaks()
    # test.test_whole_pipeline()
    exit()

    caller = test.caller
    caller._control_pileup = Pileup.from_bed_graph(
        test.graph, "control_track.bdg")
    caller._sample_pileup = Pileup.from_bed_graph(
        test.graph, "sample_track.bdg")
    cProfile.run("caller.get_score()", "profiling")
    p = pstats.Stats("profiling")
    p.sort_stats("tottime").print_stats()
    exit()
