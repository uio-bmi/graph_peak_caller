import unittest
import pytest
import numpy as np
# from graph_peak_caller.legacy.sparsepileup import SparsePileup, ValuedIndexes
# from graph_peak_caller.densepileup import DensePileup
# from graph_peak_caller.areas import BinaryContinousAreas
# from graph_peak_caller.peakscores import ScoredPeak
import offsetbasedgraph as obg
pytest.skip()


@pytest.mark.skip("Legacy")
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

        self.pileup = DensePileup.create_from_old_sparsepileup(self.pileup)

        flat_data = {node_id: ValuedIndexes(
            np.array([], dtype="int"),
            np.array([]),
            node_id*10, 10)
                     for node_id in nodes}

        self.flat_pileup = SparsePileup(self.graph)
        self.flat_pileup.data = flat_data
        self.flat_pileup = DensePileup.create_from_old_sparsepileup(self.flat_pileup)

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
            -2: ValuedIndexes(indexes, values+20, 20, 10),
            3: ValuedIndexes(indexes, values+30, 30, 10),
            -3: ValuedIndexes(indexes, values+30, 30, 10),
            4: ValuedIndexes(indexes, values+40, 40, 10),
            -4: ValuedIndexes(indexes, values+40, 40, 10),
            5: ValuedIndexes(np.array([2, 4]), np.array([52, 54]), 50, 5),
            -1: ValuedIndexes(np.array([1, 3]), np.array([16, 18]), 14, 5)
            }

        self.flat_scores = {
            2: ValuedIndexes(np.array([], dtype="int"), np.array([]), 20, 10),
            3: ValuedIndexes(np.array([], dtype="int"), np.array([]), 30, 10),
            4: ValuedIndexes(np.array([], dtype="int"), np.array([]), 40, 10),
            -2: ValuedIndexes(np.array([], dtype="int"), np.array([]), 20, 10),
            -3: ValuedIndexes(np.array([], dtype="int"), np.array([]), 30, 10),
            -4: ValuedIndexes(np.array([], dtype="int"), np.array([]), 40, 10),
            5: ValuedIndexes(np.array([], dtype="int"), np.array([]), 50, 5),
            -1: ValuedIndexes(np.array([], dtype="int"), np.array([]), 10, 5)
            }

        self.peak2 = BinaryContinousAreas(self.graph)
        self.peak2.add_full(2)
        self.peak2.add_start(-3, 7)
        self.peak2.add_full(4)
        self.peak2.add_start(5, 5)

        self.scores2 = {
            2: ValuedIndexes(indexes, values+20, 20, 10),
            -3: ValuedIndexes(np.array([1, 3, 5]), np.array([34, 36, 38]), 32, 7),
            4: ValuedIndexes(indexes, values+40, 40, 10),
            5: ValuedIndexes(np.array([2, 4]), np.array([52, 54]), 50, 5),
            }

        self.scored_peak = ScoredPeak(self.peak, self.scores)
        self.flat_scored_peak = ScoredPeak(self.peak, self.flat_scores)
        self.scored_peak2 = ScoredPeak(self.peak2, self.scores2)
        self.max_path = obg.DirectedInterval(
            5, 5, [1, 3, -4, 5],
            graph=self.graph)
        self.max_path2 = obg.DirectedInterval(
            3, 5, [3, -4, 5], graph=self.graph)

    def test_from_peak_and_pileup(self):
        scored_peak = ScoredPeak.from_peak_and_pileup(
            self.peak, self.pileup)
        self.assertEqual(scored_peak, self.scored_peak)

    def test_from_peak_and_pileup_flat(self):
        scored_peak = ScoredPeak.from_peak_and_pileup(
            self.peak, self.flat_pileup)
        #print(scored_peak)
        print(self.flat_pileup)
        self.assertEqual(scored_peak, self.flat_scored_peak)

    def test_get_max_path(self):
        max_path = self.scored_peak.get_max_path()
        self.assertEqual(max_path, self.max_path)

    def test_get_max_path_flat(self):
        max_path = self.flat_scored_peak.get_max_path()
        self.assertEqual(max_path, self.max_path)

    def test_get_max_path_two_starts(self):
        max_path = self.scored_peak2.get_max_path()
        self.assertEqual(max_path, self.max_path2)


if __name__ == "__main__":
    unittest.main()
