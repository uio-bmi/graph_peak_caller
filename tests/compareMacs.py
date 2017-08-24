from offsetbasedgraph import Graph, Block, Position, Interval
from offsetbasedgraph.interval import IntervalCollection
from graph_peak_caller.callpeaks import CallPeaks
import random


class SimpleInterval(object):
    node_id = 0

    def __init__(self, start, end, direction):
        self.start = start
        self.end = end
        self.direction = direction
        self.dir_symbol = "+" if direction > 0 else "-"

    def to_file_line(self):
        return "\t".join(
            str(v) for v in
            (self.node_id, self.start, self.end, self.dir_symbol))

    @classmethod
    def from_file_line(cls, line):
        _, start, end, dir_symbol = line.split("\t")
        direction = +1 if dir_symbol == "+" else -1
        return cls(int(start), int(end), direction)


class MACSTests(object):
    def __init__(self, node_size, n_nodes, n_intervals, read_length=50):
        self.node_size = node_size
        self.n_nodes = n_nodes
        self.n_intervals = n_intervals
        self.read_length = read_length
        self.genome_size = node_size*n_nodes

    def run(self):
        self.create_linear_graph()
        self.create_intervals()
        self.write_intervals()

        caller = CallPeaks("lin_graph", "graph_intervals")
        caller.run()

    def write_intervals(self):
        f = open("lin_intervals.bed", "w")
        f.writelines(interval.to_file_line() for
                     interval in self.linear_intervals)
        f.close()
        graph_intervals = IntervalCollection(self.graph_intervals)
        graph_intervals.to_file("graph_intervals")

    def create_linear_graph(self):
        nodes = {i: Block(self.node_size) for i in range(self.n_nodes)}
        adj_list = {i: [i+1] for i in range(self.n_nodes-1)}
        self.graph = Graph(nodes, adj_list)
        self.graph.to_file("lin_graph")

    def _get_graph_interval(self, start, end, direction):
        start_rp = start//self.node_size
        end_rp = (end+1)//self.node_size
        start_pos = Position(
            start_rp,
            start % self.node_size)
        end_pos = Position(
            end_rp,
            (end % self.node_size) + 1)
        region_paths = list(range(start_rp, end_rp+1))
        return Interval(start_pos, end_pos, region_paths, direction=direction)

    def create_intervals(self):
        self.linear_intervals = []
        self.graph_intervals = []
        for _ in range(self.n_intervals):
            direction = random.choice((-1, 1))
            start = random.randint(0, self.genome_size-1)
            end = random.randint(start+1, self.genome_size)
            self.linear_intervals.append(Interval(start, end, [0]))
            self.graph_intervals.append(
                self._get_graph_interval(start, end, direction))


if __name__ == "__main__":
    test = MACSTests(100, 100, 100)
    test.run()