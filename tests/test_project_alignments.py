from graph_peak_caller.command_line_interface import project_alignments
from offsetbasedgraph import Graph, Block, NumpyIndexedInterval, Interval

def test_simple():
    graph = Graph({1: Block(10), 2: Block(5), 3: Block(10), 4: Block(5)}, {1: [2, 3], 2: [4], 3: [4]})
    graph.convert_to_numpy_backend()
    linear_path = NumpyIndexedInterval.from_interval(Interval(0, 10, [1, 2, 4], graph))
    alignments = [Interval(5, 5, [1, 3], graph),
                  Interval(5, 5, [3, 4], graph)]
    projected = project_alignments(alignments, linear_path)
    projected = list(projected)
    assert projected[0] == (5, 15, "+")
    assert projected[1] == (15, 25, "+")


def test_advanced():

    graph = Graph({1: Block(10), 2: Block(5), 3: Block(10), 4: Block(5)}, {1: [2, 3], 2: [4], 3: [4]})
    graph.convert_to_numpy_backend()
    linear_path = NumpyIndexedInterval.from_interval(Interval(0, 10, [1, 2, 4], graph))
    alignments = [Interval(5, 7, [1, 3], graph)]
    projected = project_alignments(alignments, linear_path)
    projected = list(projected)
    assert projected[0] == (5, 17, "+")


def test_reverse():
    graph = Graph({1: Block(10), 2: Block(5), 3: Block(10), 4: Block(5)}, {1: [2, 3], 2: [4], 3: [4]})
    graph.convert_to_numpy_backend()
    linear_path = NumpyIndexedInterval.from_interval(Interval(0, 10, [1, 2, 4], graph))
    alignments = [Interval(4, 5, [-3, -1], graph)]
    projected = project_alignments(alignments, linear_path)
    projected = list(projected)
    assert projected[0] == (5, 16, "-")

