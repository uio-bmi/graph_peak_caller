from .pileup import Pileup
from .shifter import Shifter
from offsetbasedgraph.interval import IntervalCollection


class ControlTrack(object):
    def __init__(self, graph, file_name, fragment_length, extensions):
        self.graph = graph
        self.file_name = file_name
        self.fragment_length = fragment_length
        self.extensions = extensions

    def _get_pileup(self, extension):
        """TODO: read obg_alignments directly"""
        alignments = IntervalCollection.create_generator_from_file(
            self.file_name)
        shifter = Shifter(self.graph, extension)
        areas_generator = (shifter.extend_interval(alignment, 0) for alignment
                           in alignments)
        pileup = Pileup(self.graph)
        for areas in areas_generator:
            pileup.add_areas(areas)
        pileup.scale(self.fragment_length/(extension*2))
        return pileup

    def generate_background_tracks(self):
        extensions = [ext//2 for ext in self.extensions]
        return (self._get_pileup(extension) for extension in extensions)

    def generate_background_track(self):
        pileup = Pileup(self.graph, [], 0)
        genome_lambda = self.n_reads * self.fragment_length/self.genome_size
        pileup.init_value(genome_lambda)
        pileups = (self._get_pileup(extension) for extension in self.extensions)
        (pileup.update_max(ext_pileup) for ext_pileup in pileup)
        return pileup

    def create_background_tracks(self):
        for extension in self.extensions:
            shift = extension//2
            
            alignments = vg.AlignmentCollection.create_generator_from_file(self.file_name)
            obg_alignments = (alignment.path.to_obg(self.graph) for alignment in alignments)
            pileup = Pileup(self.graph, [])
            shifter = Shifter(self.graph, [], shift)
            areas_list = (shifter.extend_interval(interval, 0) for interval in obg_alignments)
            (pileup.add_areas(areas) for areas in areas_list)
            pileup.scale(extension/self.d)
            pileup.to_bed_graph(self._get_file_name(extension))
            pileup.clear()

    def _get_file_name(self, extension):
        self.file_name.replace(".json", "%sbg.bdg" % extension)

    def merge_background_tracks(self):
        pileup = Pileup(self.graph, [], 0)
        genome_lambda = control_info.n_reads *control_info * fragment_length/control_info.genome_size
        pileup.init_value(genome_lambda)
        for extension in self.extensions:
            extension_pileup = Pileup.from_file(
                self._get_file_name(extension))
            pileup.update_max(extension_pileup)

        return pileup
