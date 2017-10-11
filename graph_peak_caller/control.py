from collections import defaultdict
import numpy as np

from .pileup import Pileup
from .sparsepileup import SparsePileup
from .extender import Extender
from offsetbasedgraph.interval import IntervalCollection


class ControlTrack(object):

    def __init__(self, graph, intervals, fragment_length, extensions):
        self.graph = graph
        self.intervals = intervals
        self.fragment_length = fragment_length
        self.extensions = extensions
        self.background_pileups = []

    def _get_pileups(self, extensions):
        extenders = [Extender(self.graph, extension)
                     for extension in extensions]
        areas_generator = ((
            [extender.extend_interval(alignment, 0)
             for extender in extenders]
            for alignment in self.intervals))

        pileups = [Pileup(self.graph) for _ in extensions]

        areas_lists = [[] for _ in extensions]
        count = 0
        starts_dict_list = [defaultdict(list) for extender in extenders]
        ends_dict_list = [defaultdict(list) for extender in extenders]
        for areas in areas_generator:
            if count % 1000 == 0:
                print("#", count)
            count += 1
            for i in range(0, len(areas)):
                for node_id in areas[i].areas:
                    starts_dict_list[i][node_id].extend(areas[i].get_starts(node_id))
                    ends_dict_list[i][node_id].extend(areas[i].get_ends(node_id))
                    # areas_lists[i].append(areas[i])
                # pileups[i].add_areas(areas[i])
        # pileups = [SparsePileup.from_areas_collection(self.graph, areas)
        starts_dict_list = [{node_id: np.array(starts) for node_id, starts in starts_dict.items()}
                            for starts_dict in starts_dict_list]

        ends_dict_list = [{node_id: np.array(ends) for node_id, ends in ends_dict.items()}
                            for ends_dict in ends_dict_list]
        
        pileups = [SparsePileup.from_starts_and_ends(self.graph, starts_dict, ends_dict)
                   for starts_dict, ends_dict in zip(starts_dict_list, ends_dict_list)]

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
