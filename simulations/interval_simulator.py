import offsetbasedgraph as obg
from pyvg.sequences import SequenceRetriever
import random


class IntervalSimulator(object):
    def __init__(self, graph, interval_length):
        self._graph = graph
        self._interval_length = interval_length

    def _generate_start_positon(self):
        rp = random.choice(list(self._graph.blocks.keys()))
        offset = random.randint(0, self._graph.node_size(rp))
        return obg.Position(rp, offset)

    def generate_interval(self):
        start_pos = self._generate_start_positon()
        return self._generate_interval_from_start_pos(start_pos)

    def _generate_interval_from_start_pos(self, start_pos):
        cur_node = start_pos.region_path_id
        start_offset = start_pos.offset
        l = self._interval_length-(self._graph.node_size(cur_node)-start_offset)
        rps = [cur_node]
        while l > 0:
            cur_node = random.choice(self._graph.adj_list[cur_node])
            rps.append(cur_node)
            l -= self._graph.node_size(cur_node)

        end_offset = self._graph.node_size(cur_node) + l
        return obg.Interval(start_pos.offset, end_offset, rps)


if __name__ == "__main__":
    graph = obg.GraphWithReversals.from_file("../tests/graph.obg")
    sim = IntervalSimulator(graph, 36)
    intervals = [sim.generate_interval() for _ in range(100)]
    obg.IntervalCollection(intervals).to_file("simulated_intervals.py")
    retriever = SequenceRetriever.from_vg_graph("../tests/haplo1kg50-mhc.vg")
    sequences = [retriever.get_interval_sequence(i) for i in intervals]
    with open("simulated_sequences.fq", "w") as f:
        for i, seq in enumerate(sequences):
            f.write("@sim" + str(i) + "\n")
            f.write(seq + "\n")
            f.write("+\n")
            f.write("~"*36 + "\n")
    
