from collections import defaultdict
import numpy as np

from .pileup import Pileup
from .sparsepileup import SparsePileup
from .extender import Extender
from .areas import ValuedAreas
from offsetbasedgraph.interval import IntervalCollection


class ControlTrack(object):

    def __init__(self, graph, intervals, fragment_length, extensions):
        self.graph = graph
        self.intervals = intervals
        self.fragment_length = fragment_length
        self.extensions = extensions
        self.background_pileups = []

    def get_control_track(self, info):
        base_value = info.n_control_reads*info.fragment_length/info.genome_size
        tracks = self.generate_background_tracks()
        self.scale_pileups(tracks, base_value)
        return self.combine_backgrounds(tracks, base_value)

    def _get_pileups(self, extensions):
        extenders = [Extender(self.graph, extension)
                     for extension in extensions]
        valued_areas_list = [ValuedAreas(self.graph) for
                             _ in extensions]
        count = 0
        for alignment in self.intervals:
            if count % 1000 == 0:
                print("#", count)
            count += 1
            for extender, valued_areas in zip(extenders, valued_areas_list):
                valued_areas.add_binary_areas(
                    extender.extend_interval(alignment, 0))

        pileups = [SparsePileup.from_valued_areas(self.graph, valued_areas) for
                   valued_area in valued_areas_list]

        for pileup, extension in zip(pileups, extensions):
            pileup.scale(self.fragment_length/(extension*2))

        return pileups

    def scale_pileups(self, pileups, base_value):
        for pileup in pileups:
            pileup.scale(pileup.mean()/base_value)

    def generate_background_tracks(self):
        extensions = [ext//2 for ext in self.extensions]
        return self._get_pileups(extensions)

    def combine_backgrounds(self, background_pileups, base_value):
        max_pileup = background_pileups[0]
        for pileup in background_pileups[1:]:
            max_pileup.update_max(pileup)

        max_pileup.update_max_value(base_value)

        return max_pileup
