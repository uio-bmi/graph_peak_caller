import unittest
import numpy as np
from graph_peak_caller.sparsepileup import SparsePileup, ValuedIndexes
from graph_peak_caller.areas import BinaryContinousAreas
from graph_peak_caller.peakscores import ScoredPeak
import offsetbasedgraph as obg


class TestScoredPeak(unittest.TestCase):
    def setUp(self):
        nodes = {i+1: obg.Block(10) for i in range(5)}
        edges = {1: [2, 3], 2: [4], 3: [-4], 4: [5], -4: [5]}
        self.graph = obg.GraphWithReversals(nodes, edges)

        data = {node_id: ValuedIndexes(
            np.array([2*i for i in range(1, 5)]),
            np.array([2*i+10*node_id for i in range(1, 5)]),
            node_id*10, 10)
                for node_id in nodes}
        self.pileup = SparsePileup(self.graph)
        self.pileup.data = data

        self.peak = BinaryContinousAreas(self.graph)
        self.peak.add_full(2)
        self.peak.add_full(4)
        self.peak.add_full(3)
        self.peak.add_start(5, 5)
        self.peak.add_start(-1, 5)

        indexes = np.array([2, 4, 6, 8])
        values = np.array([2, 4, 6, 8])
        self.scores = {
            2: ValuedIndexes(indexes, values+20, 20, 10),
            3: ValuedIndexes(indexes, values+30, 30, 10),
            4: ValuedIndexes(indexes, values+40, 40, 10),
            5: ValuedIndexes(np.array([2, 4]), np.array([52, 54]), 50, 5),
            -1: ValuedIndexes(np.array([1, 3]), np.array([16, 18]), 14, 5)
            }

        self.scored_peak = ScoredPeak(self.peak, self.scores)
        self.max_path = obg.DirectedInterval(
            5, 5, [1, 3, -4, 5],
            graph=self.graph)

    def test_from_peak_and_pileup(self):
        scored_peak = ScoredPeak.from_peak_and_pileup(
            self.peak, self.pileup)
        self.assertEqual(scored_peak, self.scored_peak)

    def test_get_max_path(self):
        max_path = self.scored_peak.get_max_path()
        self.assertEqual(max_path, self.max_path)

if __name__ == "__main__":
    unittest.main()
