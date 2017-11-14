import json
import offsetbasedgraph as obg
from pybedtools import BedTool


class Peak(obg.DirectedInterval):
    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None, direction=None, score=0):
        super().__init__(start_position, end_position,
                         region_paths, graph, direction)
        self.score = score

    def set_score(self, score):
        self.score = score

    def to_file_line(self):
        object = {"start": int(self.start_position.offset),
                  "end": int(self.end_position.offset),
                  "region_paths": self.region_paths,
                  "direction": self.direction,
                  "average_q_value": self.score
                  }
        return json.dumps(object)

    @classmethod
    def from_file_line(cls, line, graph=None):
        object = json.loads(line)
        return cls(object["start"], object["end"], object["region_paths"],
                   direction=object["direction"], graph=graph)


class PeakCollection(obg.IntervalCollection):
    interval_class = Peak

    @classmethod
    def _is_in_graph(cls, peak, start_offset, end_offset):
        if peak.chrom != "chr6":
            return False
        if (peak.start < start_offset or peak.end > end_offset):
            return False
        return True

    @classmethod
    def create_from_linear_intervals_in_bed_file(
            cls, ob_graph, linear_path_interval, bed_file_name,
            graph_start_offset, graph_end_offset):
        peaks = BedTool(bed_file_name)
        intervals_on_graph = []
        i = 0
        for peak in peaks:
            start = peak.start - graph_start_offset
            end = peak.end - graph_start_offset
            if not cls._is_in_graph(peak, graph_start_offset, graph_end_offset):
                continue
            if i % 100 == 0:
                print("Interval %i" % (i))
            i += 1
            linear_interval = linear_path_interval.get_subinterval(start, end)
            linear_interval.graph = ob_graph
            intervals_on_graph.append(linear_interval)

        return cls(intervals_on_graph)

    def contains_interval(self, interval):
        for i in self.intervals:
            if i == interval:
                return True

        return False

    def get_similar_intervals(self, interval, allowed_mismatches):
        similar = []
        for i in self.intervals:
            if i.is_approx_equal(interval, allowed_mismatches):
                similar.append(i)

        return similar

    def get_identical_intervals(self, other_peak_collection):
        identical_intervals = []
        for interval in self.intervals:
            if other_peak_collection.contains_interval(interval):
                identical_intervals.append(interval)

        return identical_intervals

    def get_peaks_not_in_other_collection(self, other_collection,
                                          allowed_mismatches=0):
        out = []
        for peak in self:
            similar = other_collection.get_similar_intervals(
                peak, allowed_mismatches)
            if len(similar) == 0:
                out.append(peak)

        return out

    def get_overlapping_intervals(self, interval, minimum_overlap=1):
        overlapping = []
        for i in self.intervals:
            if i.overlaps(interval, minimum_overlap=minimum_overlap):
                overlapping.append(i)
        return overlapping
