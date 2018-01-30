from offsetbasedgraph import IntervalCollection
from collections import defaultdict
from offsetbasedgraph import IndexedInterval
import logging

class HaploTyper:
    """
    Simple naive haplotyper
    """
    def __init__(self, graph, intervals):
        assert isinstance(intervals, IntervalCollection)
        self.graph = graph
        self.intervals = intervals

        self.node_counts = defaultdict(int)

    def build(self):
        logging.info("Building haplotyper")
        i = 0
        for interval in self.intervals.intervals:
            if i % 50000 == 0:
                logging.info("%d intervals processed" % i)
            i += 1

            for rp in interval.region_paths:
                self.node_counts[rp] += 1

    def get_maximum_interval_through_graph(self):
        logging.info("Getting first blocks")
        graph = self.graph
        start = graph.get_first_blocks()
        logging.info("First blocks found")
        assert len(start) == 1, "Only works when graph has one start node"
        nodes = []
        current_block = start[0]
        i = 0
        while True:
            if i % 1000000 == 0:
                logging.info("Processing node %d" % i)
            i += 1
            nodes.append(current_block)
            next_blocks = graph.adj_list[current_block]
            if len(next_blocks) < 1:
                break

            next = None
            max_pileup_value = -1
            for potential_next in next_blocks:
                pileup_value = self.node_counts[potential_next]
                if pileup_value > max_pileup_value:
                    next = potential_next
                    max_pileup_value = pileup_value

            current_block = next

        return IndexedInterval(0, graph.node_size(nodes[-1]), nodes, graph)
