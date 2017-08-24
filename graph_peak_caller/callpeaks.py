from offsetbasedgraph import IntervalCollection
import offsetbasedgraph
import numpy as np
import pyvg as vg
from graph_peak_caller import Shifter
from graph_peak_caller import get_shift_size_on_offset_based_graph
from .control import ControlTrack
from .bdgcmp import *


class CallPeaks(object):
    def __init__(self, graph_file_name, sample_file_name,
                 control_file_name=None):
        self.graph_file_name = graph_file_name
        self.sample_file_name = sample_file_name
        self.control_file_name = control_file_name if control_file_name is not None else sample_file_name
        self.has_control = control_file_name is not None
        self._p_value_track = "p_value_track"
        self.read_length = 42

    def run(self):
        self.create_graph()
        self.find_info()
        self.determine_shift()
        self.sample_file_name = self.remove_alignments_not_in_graph(self.sample_file_name)
        if self.control_file_name is not None:
            self.control_file_name = self.remove_alignments_not_in_graph(self.control_file_name)

        self.sample_file_name = self.filter_duplicates(self.sample_file_name)
        if self.control_file_name is not None:
            self.control_file_name = self.filter_duplicates(self.control_file_name)

        self.create_control()
        self.create_sample_pileup()
        self.get_p_values()
        self.call_peaks()

    def remove_alignments_not_in_graph(self, alignment_file_name):
        interval_collection = IntervalCollection.create_generator_from_file(alignment_file_name)
        filtered_file_name = alignment_file_name + "_filtered"
        filtered_file = open(filtered_file_name, "w")
        for interval in self._get_intervals_in_ob_graph(interval_collection):
            filtered_file.writelines(["%s\n" % interval.to_file_line()])
        filtered_file.close()
        print("Alignments without duplicates written to %s" % filtered_file_name)
        return filtered_file_name

    def filter_duplicates(self, alignment_file_name):
        interval_collection = IntervalCollection.create_generator_from_file(alignment_file_name)
        filtered_file_name = alignment_file_name + "_filtered_duplicates"
        filtered_file = open(filtered_file_name, "w")

        interval_hashes = {}
        for interval in interval_collection:
            hash = interval.hash()
            if hash in interval_hashes:
                #print("Duplicate found")
                continue

            interval_hashes[hash] = True
            filtered_file.writelines(["%s\n" % interval.to_file_line()])
        filtered_file.close()
        print("Filtered alignments written to %s" % filtered_file_name)
        return filtered_file_name

    def _get_intervals_in_ob_graph(self, intervals):
        # Returns only those intervals that exist in graph
        for interval in intervals:
             if interval.region_paths[0] in self.ob_graph.blocks:
                 yield interval

    def find_info(self):
        genome_size = 0
        #lines = (line["node"] for line in self.graph_file_name.readlines() if "node" in line)
        #sizes = (sum(Node.from_json(json_obj).n_basepairs for json_obj in line) for line in lines)
        sizes = (block.length() for block in self.ob_graph.blocks.values())

        self.genome_size = sum(sizes)
        self.n_reads = sum(1 for line in open(self.control_file_name))

    def determine_shift(self):
        try:
            self.shift = get_shift_size_on_offset_based_graph(self.ob_graph, self.sample_file_name)
        except RuntimeError:
            print("WARNING: To litle data to compute shift. Setting default shift")
            self.shift = 125

    def create_graph(self):
        #self.vg_graph = vg.Graph.create_from_file(self.graph_file_name)
        self.ob_graph = offsetbasedgraph.Graph.from_file(self.graph_file_name)

    def create_control(self):
        extensions = [self.shift, 1000, 5000, 10000] if self.has_control else [5000, 10000]
        control_track = ControlTrack(
            self.ob_graph, self.control_file_name,
            self.shift, extensions)
        background_value = self.n_reads*self.shift/self.genome_size
        self._control_track = create_background_pileup_as_max_from_pileups(
            self.ob_graph, control_track.generate_background_tracks(),
            background_value, "control_track.bdg")

    def get_p_values(self):
        self.p_values = get_p_value_track(self.ob_graph, self._control_track, self._sample_track, self._p_value_track)

    def call_peaks(self, cutoff=0.05):
        self.p_values.threshold(-np.log10(cutoff))
        self.p_values.fill_small_wholes(self.read_length)
        self.final_track = self.p_values
        self.final_track.to_bed_graph("final_track")

    def create_sample_pileup(self):
        alignments = IntervalCollection.create_generator_from_file(
            self.sample_file_name)

        shifter = Shifter(self.ob_graph, self.shift)
        areas_list = (shifter.extend_interval(interval)
                      for interval in alignments)
        pileup = Pileup(self.ob_graph)
        for areas in areas_list:
            pileup.add_areas(areas)
        self._sample_track = "sample_track.bdg"
        pileup.to_bed_graph(self._sample_track)

    def _write_vg_alignments_as_intervals_to_bed_file(self):
        pass


if __name__ == "__main__":

    chromosome = "chr2R"
    vg_graph = vg.Graph.create_from_file("dm_test_data/x_%s.json" % chromosome, 30000, chromosome)
    ofbg = vg_graph.get_offset_based_graph()
    interval_file = vg.util.vg_mapping_file_to_interval_file("intervals_reads3_chr2R", vg_graph, "dm_test_data/reads3.json", ofbg)
    ofbg.to_file("graph.tmp")

    caller = CallPeaks("graph.tmp", interval_file)
    caller.create_graph()
    caller.find_info()
    caller.determine_shift()
    caller.sample_file_name = caller.remove_alignments_not_in_graph(caller.sample_file_name)
    caller.sample_file_name = caller.filter_duplicates(caller.sample_file_name)
