from .pileup import Pileup
from .sparsepileup import SparsePileup
from .shifter import Shifter
from offsetbasedgraph.interval import IntervalCollection


class ControlTrack(object):

    def __init__(self, graph, intervals, fragment_length, extensions):
        self.graph = graph
        self.intervals = intervals
        self.fragment_length = fragment_length
        self.extensions = extensions
        self.background_pileups = []

    def _get_pileups(self, extensions):
        shifters = [Shifter(self.graph, extension) for extension in extensions]
        areas_generator = ((
            [shifter.extend_interval_fast(alignment, 0)
             for shifter in shifters]
            for alignment in self.intervals))

        pileups = [Pileup(self.graph) for _ in extensions]

        areas_lists = [[] for _ in extensions]
        for areas in areas_generator:
            for i in range(0, len(areas)):
                areas_lists[i].append(areas[i])
                # pileups[i].add_areas(areas[i])
        pileups = [SparsePileup.from_areas_collection(self.graph, areas)
                   for areas in areas_lists]

        for i in range(0, len(pileups)):
            pileups[i].scale(self.fragment_length/(extensions[i]*2))

        return pileups

    def generate_background_tracks(self):
        extensions = [ext//2 for ext in self.extensions]
        return self._get_pileups(extensions)

    def _combine_backgrounds(self, background_pileups, base_value):
        pileup = Pileup(self.graph)
        pileup.init_value(base_value)
        for new_pileup in background_pileups:
            pileup.update_max(new_pileup)
        return pileup

    def combine_backgrounds(self, background_pileups, base_value):
        max_pileup = background_pileups[0]
        for pileup in background_pileups[1:]:
            max_pileup.update_max(pileup)

        max_pileup.update_max_value(base_value)

        return max_pileup
