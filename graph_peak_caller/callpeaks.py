from offsetbasedgraph import IntervalCollection
import pyvg as vg
from graph_peak_caller import Shifter
from graph_peak_caller import get_shift_size_on_offset_based_graph

class CallPeaks(object):
    def __init__(self, graph_file_name, sample_file_name, control_file_name=None):
        self.graph_file_name = graph_file_name
        self.sample_file_name = sample_file_name
        self.control_file_name = control_file_name if control_file_name is not None else sample_file_name
        self.has_control = control_file_name is not None

    def run(self):
        self.find_info()
        self.determine_shift()
        self.create_graphs()

        self.sample_file_name = self.remove_alignments_not_in_graph(self.sample_file_name)
        if self.control_file_name is not None:
            self.control_file_name = self.remove_alignments_not_in_graph(self.control_file_name)

        self.sample_file_name = self.filter_duplicates(self.sample_file_name)
        if self.control_file_name is not None:
            self.control_file_name = self.filter_duplicates(self.control_file_name)


        self.create_control()
        self.create_sample_pileup()
        self.get_p_values()
        self.find_peaks()

    def remove_alignments_not_in_graph(self, alignment_file_name):
        interval_collection = IntervalCollection.create_generator_from_file(alignment_file_name)
        filtered_file_name = alignment_file_name + "_filtered"
        filtered_file = open(filtered_file_name, "w")
        for interval in self._get_intervals_in_ob_graph(interval_collection):
            filtered_file.writelines(["%s\n" % interval.to_file_line()])
        filtered_file.close()
        return filtered_file_name

    def filter_duplicates(self, alignment_file_name):
        interval_collection = IntervalCollection.create_generator_from_file(alignment_file_name)
        filtered_file_name = alignment_file_name + "_filtered_duplicates"
        filtered_file = open(filtered_file_name, "w")

        interval_hashes = {}
        for interval in interval_collection:
            hash = interval.hash()
            if hash in interval_hashes:
                print("Duplicate found")
                continue

            interval_hashes[hash] = True
            filtered_file.writelines(["%s\n" % interval.to_file_line()])
        filtered_file.close()
        return filtered_file_name

    def _get_intervals_in_ob_graph(self, intervals):
        # Returns only those intervals that exist in vg_graph
        for interval in intervals:
             if interval.region_paths[0] in self.vg_graph.blocks:
                 yield interval

    def find_info(self):
        genome_size = 0
        lines = (line["node"] for line in self.graph_file_name.readlines() if "node" in line)
        sizes = (sum(Node.from_json(json_obj).n_basepairs for json_obj in line) for line in lines)

        self.genome_size = sum(sizes)
        self.n_reads = sum(1 for line in open(self.control_file_name))

    def determine_shift(self):
        self.shift = get_shift_size_on_offset_based_graph(self.vg_graph, self.sample_file_name, self.chromosome, self.linear_genome_size)

    def create_graphs(self):
        self.vg_graph = vg.Graph.create_from_file(self.graph_file_name)
        self.ob_graph = self.vg_graph.get_offset_based_graph()

    def create_control(self):
        bg_track = BackgroundTrack(self.ob_graph, control_file_name, self.d, 1000, 10000)
        f = open()
        jsons = (json.loads(line) for line in f.readlines())
        alignments = [vg.Alignment.from_json(json_dict) for json_dict in jsons]

    def create_sample_pileup(self):
        alignments = vg.AlignmentCollection.create_generator_from_file(
            self.control_file_name)
        obg_alignments = (alignment.path.to_obg(ob_graph) for alignment in alignments)
        shifter = Shifter(self.ob_graph, obg_alignments, self.shift)
        areas_list = (shifter.extend_interval(interval) for interval in obg_alignments)
        pileup = Pileup(self.ob_graph, [])
        (pileup.add_areas() for area in areas)

    def _write_vg_alignments_as_intervals_to_bed_file(self):
        pass

