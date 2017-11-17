import unittest
from offsetbasedgraph import GraphWithReversals, Block, Interval, DirectedInterval
from graph_peak_caller.callpeaks import ExperimentInfo, CallPeaks
# Testing of all steps before peak calling in callpeaks

class TestPreCallpeaks(unittest.TestCase):

    def setUp(self):

        self.fragment_length = 5
        self.read_length = 2

        self.split_graph = GraphWithReversals(
            {
                1: Block(10),
                2: Block(10)
            }
        )

        pass

    def test_simple_without_control(self):



