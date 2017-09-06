from .pileup import Pileup
from .shifter import Shifter
from offsetbasedgraph.interval import IntervalCollection


class ControlTrack(object):
    def __init__(self, graph, intervals, fragment_length, extensions):
        self.graph = graph

        self.intervals = intervals
        assert isinstance(intervals, str), "Intervals must be a file name"
        #if isinstance(intervals, str):
        #    self.intervals = IntervalCollection.create_generator_from_file(intervals)
        self.fragment_length = fragment_length
        self.extensions = extensions

    def _get_pileup(self, extension):
        """TODO: read obg_alignments directly"""
        alignments = IntervalCollection.from_file(self.intervals)
        shifter = Shifter(self.graph, extension)
        areas_generator = (shifter.extend_interval_fast(alignment, 0) for alignment
                           in alignments)
        pileup = Pileup(self.graph)
        for areas in areas_generator:
            pileup.add_areas(areas)
        pileup.scale(self.fragment_length/(extension*2))
        return pileup

    def generate_background_tracks(self):
        extensions = [ext//2 for ext in self.extensions]
        return (self._get_pileup(extension) for extension in extensions)

    def combine_backgrounds(self, backgrounds, base_value):
        pileup = Pileup(self.graph)
        pileup.init_value(base_value)
        for new_pileup in backgrounds:
            pileup.update_max(new_pileup)
        return pileup

    def _get_file_name(self, extension):
        return "intervals_%sbg.bdg" % extension
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
