import unittest
from graph_peak_caller.subgraphcollection import SubgraphCollection, ConnectedAreas
from graph_peak_caller.extender import Areas
from offsetbasedgraph import GraphWithReversals as Graph, Block, \
    Interval, Position, GraphWithReversals
import numpy as np
#from graph_peak_caller.sparsepileupv2 import SparsePileup, SparsePileupData
from graph_peak_caller.densepileup import DensePileup
from graph_peak_caller.subgraphcollection import \
    SubgraphCollectionPartiallyOrderedGraphBuilder, SubgraphCollectionPartiallyOrderedGraph
from graph_peak_caller.peakscores import ScoredPeak

from graph_peak_caller.areas import BinaryContinousAreas

class SubgraphsAndMaxPathsOnPartiallyOrderedGraphs(unittest.TestCase):

    def setUp(self):
        self.linear_graph = Graph({i: Block(5) for i in range(1, 4)},
                                  {i: [i+1] for i in range(1, 3)})

        self.scores = DensePileup.from_intervals(self.linear_graph,
                                                  [Interval(0, 5, [i]) for i in range(1, 4)])

        self.graph = Graph(
            {i: Block(5) for i in range(1, 4)},
            {
                1: [3],
                2: [3],
                3: [4]
            }
        )

    def test_simple(self):

        intervals = [Interval(1, 2, [1, 2]), Interval(1, 4, [3])]
        pileup = DensePileup.from_intervals(self.linear_graph,
                                             intervals)

        subgraphs = SubgraphCollectionPartiallyOrderedGraph.create_from_pileup(
            self.linear_graph, pileup)

        scored_peaks = (ScoredPeak.from_peak_and_pileup(peak, self.scores)
                        for peak in subgraphs)
        max_paths = [peak.get_max_path() for peak in scored_peaks]
        self.assertTrue(all(interval in max_paths for interval in intervals))


    def test_simple2(self):
        intervals = [Interval(1, 5, [1]), Interval(1, 2, [2, 3])]
        pileup = DensePileup.from_intervals(self.graph,
                                             intervals)

        subgraphs = SubgraphCollectionPartiallyOrderedGraph.create_from_pileup(
            self.graph, pileup)
        #print(subgraphs)

        scored_peaks = (ScoredPeak.from_peak_and_pileup(peak, self.scores)
                        for peak in subgraphs)
        max_paths = [peak.get_max_path() for peak in scored_peaks]
        print(max_paths)
        self.assertTrue(Interval(1, 2, [1, 3]) in max_paths or
                        Interval(1, 2, [2, 3]) in max_paths)

    def test_simple3(self):
        graph = Graph({i: Block(5) for i in range(1, 6)},
                      {
                          1: [3],
                          2: [3],
                          3: [4, 5]
                      })
        scores = DensePileup.from_intervals(graph,
                                            [Interval(0, 5, [i]) for i in range(1, 6)])
        intervals = [Interval(0, 5, [1]),
                     Interval(0, 5, [3]),
                     Interval(0, 5, [4]),
                     Interval(0, 3, [5])]
        pileup = DensePileup.from_intervals(graph,
                                             intervals)
        subgraphs = SubgraphCollectionPartiallyOrderedGraph.create_from_pileup(
            graph, pileup)
        scored_peaks = (ScoredPeak.from_peak_and_pileup(peak, scores)
                        for peak in subgraphs)
        max_paths = [peak.get_max_path() for peak in scored_peaks]

        self.assertTrue(
            Interval(0, 5, [1, 3, 4]) in max_paths or
            Interval(0, 3, [1, 3, 5]) in max_paths
        )

if __name__ == "__main__":
    unittest.main()