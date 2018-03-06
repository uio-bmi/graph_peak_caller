import json
import offsetbasedgraph as obg
import logging
from collections import defaultdict
import numpy as np
from .nongraphpeaks import NonGraphPeakCollection, NonGraphPeak
from offsetbasedgraph import DirectedInterval


class Peak(obg.DirectedInterval):
    def __init__(self, start_position, end_position,
                 region_paths=None, graph=None, direction=None, score=0, unique_id=None):
        super().__init__(start_position, end_position,
                         region_paths, graph, direction)
        self.score = score
        self.unique_id = unique_id
        self.info = (0, 0)

    def __str__(self):
        return super().__str__() + " [%s]" % self.score

    @classmethod
    def from_interval_and_score(cls, interval, score):
        return cls(interval.start_position, interval.end_position,
                   interval.region_paths, interval.graph, interval.direction,
                   score)

    def set_score(self, score):
        assert isinstance(score, float), "Score %s is invalid, type %s" % (score, type(score))
        self.score = score

    def to_file_line(self):
        object = {"start": int(self.start_position.offset),
                  "end": int(self.end_position.offset),
                  "region_paths": [int(r) for r in self.region_paths],
                  "direction": int(self.direction),
                  "average_q_value": float(self.score),
                  "unique_id": str(self.unique_id),
                  "is_ambigous": int(self.info[1]),
                  "is_diff": int(self.info[0])
                  }
        try:
            d = json.dumps(object)
        except:
            for k, v in object.items():
                print(k, v, type(v))
            raise
        return d

    @classmethod
    def from_file_line(cls, line, graph=None):
        object = json.loads(line)
        unique_id = None
        if "unique_id" in object:
            unique_id = object["unique_id"]
        obj = cls(object["start"], object["end"], object["region_paths"],
                  direction=object["direction"], graph=graph,
                  score=object["average_q_value"],
                  unique_id=unique_id)
        obj.info = (object["is_diff"], object["is_ambigous"])
        return obj

    def to_approx_linear_peak(self, linear_path, chromosome):
        # Assuming graph is partially ordered DAG
        intersecting_nodes = set(self.region_paths).intersection(linear_path.nodes_in_interval())
        intersecting_nodes = sorted(list(intersecting_nodes))

        first_node = intersecting_nodes[0]
        start_offset = 0
        if first_node == self.region_paths[0]:
            start_offset = self.start_position.offset

        last_node = intersecting_nodes[-1]
        end_offset = 0
        if last_node == self.region_paths[-1]:
            end_offset = self.end_position.offset

        linear_start_pos = linear_path.get_offset_at_node(first_node) + start_offset
        linear_end_pos = linear_path.get_offset_at_node(last_node) + end_offset
        return NonGraphPeak(chromosome, linear_start_pos, linear_end_pos)


class PeakCollection(obg.IntervalCollection):
    interval_class = Peak

    @classmethod
    def _is_in_graph(cls, peak, chrom, start_offset, end_offset):
        if peak.chromosome != chrom:
            return False
        if (peak.start < start_offset or peak.end > end_offset):
            return False
        return True

    def create_node_index(self):
        # Creates an index from region path to intervals touching the rp,
        # making it fast to do approx. overlap
        index = defaultdict(list)
        for peak in self.intervals:
            for rp in peak.region_paths:
                index[abs(rp)].append(peak)

        self._index = index

    def approx_contains_part_of_interval(self, interval):
        assert hasattr(self, "_index"), "Create index first by calling create_node_index()"
        for rp in interval.region_paths:
            if len(self._index[abs(rp)]) > 0:
                return True

        return False

    def which_approx_contains_part_of_interval(self, interval):
        assert hasattr(self, "_index"), "Create index first by calling create_node_index()"
        for rp in interval.region_paths:
            containing = self._index[abs(rp)]
            if len(containing) > 0:
                return containing[0]

    @classmethod
    def create_from_nongraph_peak_collection(
            cls, ob_graph, peak_collection,
            linear_path_interval,
            graph_region=None):
        peaks = peak_collection.peaks
        intervals_on_graph = []
        i = 0
        graph_start_offset = 0 if graph_region is None else graph_region.start
        for peak in peaks:
            start = peak.start - graph_start_offset
            end = peak.end - graph_start_offset
            end = min(end, linear_path_interval.length())
            if graph_region is not None:
                if not cls._is_in_graph(peak, graph_region.chromosome,
                                        graph_region.start,
                                        graph_region.end):
                    logging.info("Filtered out peak")
                    continue
            if i % 100 == 0:
                logging.info("%d peaks processed" % (i))
            i += 1
            linear_interval = linear_path_interval.get_exact_subinterval(start, end)
            linear_interval.graph = ob_graph
            graph_peak = Peak.from_interval_and_score(
                linear_interval, peak.score)
            graph_peak.unique_id = peak.unique_id
            graph_peak.sequence = peak.sequence
            intervals_on_graph.append(graph_peak)

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

    def to_approx_linear_peaks(self, linear_path, chromosome):
        linear_peaks = []
        i = 0
        for peak in self.intervals:
            if i % 500 == 0:
                logging.info("Converting peak %d" % i)
            i += 1
            linear_peaks.append(peak.to_approx_linear_peak(linear_path, chromosome))

        return NonGraphPeakCollection(linear_peaks)


    @classmethod
    def from_fasta_file(cls, file_name, graph=None):
        f = open(file_name)
        peaks = []
        while True:
            header = f.readline()
            sequence = f.readline()
            if not sequence:
                break

            header = header.split(maxsplit=1)
            id = header[0].replace(">", "")
            interval_json = header[1]
            peak = Peak.from_file_line(interval_json)
            peak.unique_id = id
            peak.sequence = sequence
            peak.graph = graph
            peaks.append(peak)
        avg_peak_size = np.mean([p.length() for p in peaks])
        logging.info("Avg peak size: %2.f" % avg_peak_size)

        return cls(peaks)

