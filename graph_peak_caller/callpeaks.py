from offsetbasedgraph import IntervalCollection
import offsetbasedgraph
import numpy as np
import pyvg as vg
from graph_peak_caller import Shifter
from graph_peak_caller import get_shift_size_on_offset_based_graph
from .control import ControlTrack
from .bdgcmp import *


def enable_filewrite(func):
    def wrapper(*args, **kwargs):
        intervals = args[1]
        if isinstance(intervals, str):
            intervals = IntervalCollection.from_file(intervals)

        write_to_file = kwargs.pop("write_to_file", False)
        interval_list = func(args[0], intervals, **kwargs)

        if write_to_file:
            interval_collection = IntervalCollection(interval_list)
            interval_collection.to_file(write_to_file)
            return write_to_file

            with open(write_to_file, "w") as file:
                print("Wrote results to " + str(write_to_file))
                file.writelines(("%s\n" % interval.to_file_line() for interval in interval_list))

            return write_to_file
        else:
            return interval_list

    return wrapper


class ExperimentInfo(object):
    def __init__(self, genome_size, n_sample_reads,
                 n_control_reads, fragment_length, read_length):
        self.genome_size = genome_size
        self.n_sample_reads = n_sample_reads
        self.n_control_reads = n_control_reads
        self.fragment_length = fragment_length
        self.read_length = read_length

    @classmethod
    def find_info(cls, graph, sample_file_name, control_file_name=None):
        sizes = (block.length() for block in graph.blocks.values())
        genome_size = sum(sizes)
        n_sample_reads = sum(1 for line in open(control_file_name))
        n_control_reads = n_sample_reads
        if control_file_name is not None:
            n_control_reads = sum(1 for line in open(control_file_name))
        try:
            fragment_length, read_length = get_shift_size_on_offset_based_graph(
                graph, sample_file_name)
            print("Found fragment length=%d, read length=%d" % (fragment_length, read_length))
        except RuntimeError:
            print("WARNING: To liptle data to compute shift. Setting to default.")
            fragment_length = 125
            read_length = 20
        return cls(genome_size, n_sample_reads, n_control_reads, fragment_length, read_length)


class CallPeaks(object):
    def __init__(self, graph_file_name, sample_file_name,
                 control_file_name=None, experiment_info=None, verbose=False):
        self.graph_file_name = graph_file_name
        self.sample_file_name = sample_file_name
        self.has_control = control_file_name is not None
        self.control_file_name = control_file_name if self.has_control else sample_file_name
        self._p_value_track = "p_value_track"
        self._q_value_track = "q_value_track"
        self.info = experiment_info
        self.verbose = verbose
        self.sample_intervals = []
        self.control_intervals = []
        self._control_pileup = None
        self._sample_pileup = None

    def run(self, out_file="final_peaks"):
        self.create_graph()
        self.preprocess()
        if self.info is None:
            self.info = ExperimentInfo.find_info(
                self.ob_graph, self.sample_file_name, self.control_file_name)
        self.create_control()
        self.create_sample_pileup()
        self.scale_tracks()
        self.get_p_values()
        self.get_q_values()
        self.call_peaks(out_file)

    def preprocess(self):
        self.sample_intervals = self.remove_alignments_not_in_graph(self.sample_file_name)
        self.sample_file_name = self.filter_duplicates(self.sample_intervals, write_to_file=self.sample_file_name+"_filtered")

        if self.control_file_name is not None:
            self.control_intervals = self.remove_alignments_not_in_graph(self.control_file_name)
            self.control_file_name = self.filter_duplicates(self.control_intervals, write_to_file=self.control_file_name+"_filtered")

        self.info.n_sample_reads = self.count_number_of_intervals_in_file(self.sample_file_name)
        self.info.n_control_reads = self.count_number_of_intervals_in_file(self.control_file_name)

    def count_number_of_intervals_in_file(self, interval_file_name):
        n = 0
        for interval in IntervalCollection.from_file(interval_file_name):
            n += 1
        print("Number of intervals: %d" % n)
        return n

    @enable_filewrite
    def remove_alignments_not_in_graph(self, intervals):
        for interval in self._get_intervals_in_ob_graph(intervals):
            yield interval

    @enable_filewrite
    def filter_duplicates(self, intervals):

        interval_hashes = {}
        n_duplicates = 0
        n_reads_left = 0
        for interval in intervals:
            hash = interval.hash()
            if hash in interval_hashes:
                n_duplicates += 1
                continue

            interval_hashes[hash] = True

            yield interval

    def _get_intervals_in_ob_graph(self, intervals):
        # Returns only those intervals that exist in graph
        for interval in intervals:
            if interval.region_paths[0] in self.ob_graph.blocks:
                yield interval

    def scale_tracks(self, update_saved_files=False):
        print("Scaling tracks to ratio: %d / %d" % (self.info.n_sample_reads, self.info.n_control_reads))
        ratio = self.info.n_sample_reads/self.info.n_control_reads

        if self.info.n_sample_reads == self.info.n_control_reads:
            return

        if ratio > 1:
            self._sample_pileup.scale(1/ratio)
            if update_saved_files:
                self._sample_pileup.to_bed_graph(self._sample_track)
        else:
            self._control_pileup(ratio)
            if update_saved_files:
                self._control_pileup.to_bed_graph(self._control_track)

        #scale_down_tracks(ratio, self.sample_file_name, self.control_file_name)

    def find_info(self):
        genome_size = 0
        #lines = (line["node"] for line in self.graph_file_name.readlines() if "node" in line)
        #sizes = (sum(Node.from_json(json_obj).n_basepairs for json_obj in line) for line in lines)
        sizes = (block.length() for block in self.ob_graph.blocks.values())

        self.genome_size = sum(sizes)
        self.n_reads = sum(1 for line in open(self.control_file_name))

    def create_graph(self):
        self.ob_graph = offsetbasedgraph.Graph.from_file(self.graph_file_name)

    def create_control(self, save_to_file=False):
        if self.verbose:
            print("Creating control")
        extensions = [self.info.fragment_length, 2500, 5000] if self.has_control else [5000]

        control_track = ControlTrack(
            self.ob_graph, self.control_file_name,
            self.info.fragment_length, extensions)

        tracks = control_track.generate_background_tracks()
        background_value = self.info.n_control_reads*self.info.fragment_length/self.info.genome_size
        print("Background value: %.4f" % background_value)
        print("N control reads: %d"  % self.info.n_control_reads)
        print(self.info.fragment_length)
        print(self.info.genome_size)

        pileup = control_track.combine_backgrounds(tracks, background_value)

        if save_to_file:
            self._control_track = "control_track.bdg"
            pileup.to_bed_graph(self._control_track)
            print("Saved control pileup to " + self._control_track)

        self._control_pileup = pileup

    def get_q_values(self):
        get_q_values_track_from_p_values(self.p_values)

    def get_p_values(self):
        print("Get p-values")
        self.p_values = get_p_value_track_from_pileups(self.ob_graph, self._control_pileup, self._sample_pileup)
        #self.p_values.to_bed_graph(self._p_value_track)
        print(self.p_values)

    def call_peaks(self, out_file="final_peaks", cutoff=0.05):
        print("Calling peaks")
        self.p_values.threshold(-np.log10(cutoff))
        self.p_values.fill_small_wholes(self.info.read_length)
        self.p_values.remove_small_peaks(self.info.fragment_length)
        self.final_track = self.p_values
        self.final_track.to_bed_file(out_file)

    def create_sample_pileup(self, save_to_file=False):
        print("Create sample pileup")
        alignments = IntervalCollection.from_file(
            self.sample_file_name)

        shifter = Shifter(self.ob_graph, self.info.fragment_length)
        areas_list = (shifter.extend_interval_fast(interval)
                      for interval in alignments)
        pileup = Pileup(self.ob_graph)
        for areas in areas_list:
            pileup.add_areas(areas)
        self._sample_track = "sample_track.bdg"
        if save_to_file:
            pileup.to_bed_graph(self._sample_track)
            print("Saved sample pileup to " + self._sample_track)

        self._sample_pileup = pileup

    def _write_vg_alignments_as_intervals_to_bed_file(self):
        pass


if __name__ == "__main__":
    chromosome = "chr2R"
    vg_graph = vg.Graph.create_from_file("dm_test_data/x_%s.json" % chromosome, 30000, chromosome)
    ofbg = vg_graph.get_offset_based_graph()
    interval_file = vg.util.vg_mapping_file_to_interval_file("intervals_reads3_chr2R", vg_graph, "dm_test_data/reads3_small.json", ofbg)
    ofbg.to_file("graph.tmp")

    caller = CallPeaks("graph.tmp", interval_file)
    caller.create_graph()
    caller.find_info()
    caller.determine_shift()
    caller.sample_file_name = caller.remove_alignments_not_in_graph(caller.sample_file_name)
    caller.sample_file_name = caller.filter_duplicates(caller.sample_file_name)
