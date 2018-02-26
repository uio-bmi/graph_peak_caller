from pyvg.conversion import vg_json_file_to_intervals
import offsetbasedgraph as obg
from graph_peak_caller.haplotyping import HaploTyper
import logging


def graph_to_indexed_interval(graph):
    pass


class LinearFilter:
    def __init__(self, position_tuples, indexed_interval):
        self._nodes = set()
        self._nodes.update(indexed_interval.region_paths)
        self._position_tuples = position_tuples
        self._indexed_interval = indexed_interval

    def find_start_positions(self):
        start_positions = {"+": [],
                 "-": [] }

        logging.info("Mapping graph position to linear positions")
        nodes_in_linear_path = set(self._indexed_interval.region_paths)
        i = 0
        for pos in self._position_tuples:
            direction = "+"
            node = pos.region_path_id
            if node  < 0:
                direction = "-"

            if i % 10000 == 0:
                logging.info("%d reads processed" % i)
            i += 1

            if abs(node) not in nodes_in_linear_path:
                continue

            offset = self._indexed_interval.get_offset_at_position(
                pos, direction)
            assert isinstance(offset, int), "Offset is of type %s" % type(offset)
            start_positions[direction].append(offset)

        logging.info("Found in total %d positions" % len(start_positions["+"]))
        return start_positions

    @classmethod
    def from_vg_json_reads_and_graph(cls, json_file_name, graph_file_name):
        logging.info("Reading graph %s" % graph_file_name)
        graph = obg.GraphWithReversals.from_unknown_file_format(graph_file_name)

        logging.info("Getting indexed interval through graph")
        intervals =  vg_json_file_to_intervals(None, json_file_name, graph)
        haplotyper = HaploTyper(graph, obg.IntervalCollection(intervals))
        haplotyper.build()
        indexed_interval = haplotyper.get_maximum_interval_through_graph()
        #indexed_interval = graph.get_indexed_interval_through_graph()

        intervals =  vg_json_file_to_intervals(None, json_file_name, graph)
        positions = (interval.start_position for interval in intervals)

        return cls(positions, indexed_interval)

