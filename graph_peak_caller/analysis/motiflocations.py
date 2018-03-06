from collections import defaultdict


class MotifLocations:
    def __init__(self, peaks, motifs):
        self._peaks = {peak.unique_id: peak for peak in peaks}
        self._motifs = motifs
        self._locations = self._get_motif_locations()

    def _get_motif_locations(self):
        locations = []
        for peak_id, interval in self._peaks.items():
            motif_locations = self._motifs.get_entries(peak_id)
            for motif_location in motif_locations:
                locations.append(
                    interval.get_subinterval(motif_location._start,
                                             motif_location._end))
        return locations

    def get_rp_index(self):
        index = defaultdict(list)
        for location in self._locations:
            for rp in location.region_paths:
                index[rp].append(location)

    def find_matches(self, other):
        index = self.get_rp_index()
        match_dict = {}
        for location in other._locations:
            matches = []
            for rp in location.region_paths:
                matches.extend(index[rp])
            match_dict[location] = matches
        return match_dict
