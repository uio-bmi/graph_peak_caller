from offsetbasedgraph.interval import IntervalCollection
from pybedtools import BedTool


class PeakCollection(IntervalCollection):

    def __init__(self, interval_list):
        super(PeakCollection, self).__init__(interval_list)

    @classmethod
    def create_from_linear_intervals_in_bed_file(cls,
                        ob_graph, linear_path_interval, bed_file_name, graph_start_offset, graph_end_offset):
        peaks = BedTool(bed_file_name)
        intervals_on_graph = []
        i = 0
        for peak in peaks:
            start = peak.start - graph_start_offset
            end = peak.end - graph_start_offset
            if peak.chrom != "chr6" or peak.start < graph_start_offset or peak.end > graph_end_offset:
                continue
            #print("Peak  %i" % (i))
            i += 1
            #if i > 10:
            #    break
            linear_interval = linear_path_interval.get_subinterval(start, end)
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
            if other_peak_collection.contains(interval):
                identical_intervals.append()


    def get_peaks_not_in_other_collection(self, other_collection, allowed_mismatches = 0):
        out = []
        for peak in self:
            similar = other_collection.get_similar_intervals(peak, allowed_mismatches)
            if len(similar) == 0:
                out.append(peak)

        return out
