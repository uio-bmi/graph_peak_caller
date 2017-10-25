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
        for peak in peaks:
            start = peak.start - graph_start_offset
            end = peak.end - graph_start_offset
            if peak.start < graph_start_offset or peak.end > graph_end_offset:
                continue
            print(start, end)
            print("linear interval")
            linear_interval = linear_path_interval.get_subinterval(start, end)
            print("Linear interval")
            intervals_on_graph.append(linear_interval)

        return cls(intervals_on_graph)

    def contains_interval(self, interval):
        for i in self.intervals:
            if i == interval:
                return True

        return False

    def get_identical_intervals(self, other_peak_collection):
        identical_intervals = []
        for interval in self.intervals:
            if other_peak_collection.contains(interval):
                identical_intervals.append()






