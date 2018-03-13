import unittest
import pytest
pytestmark = pytest.mark.skip
# from graph_peak_caller.subgraphcollection import SubgraphCollection, ConnectedAreas
from offsetbasedgraph import GraphWithReversals as Graph, Block, Interval
from graph_peak_caller.areas import BinaryContinousAreas, Areas
from graph_peak_caller.densepileup import DensePileup
# from graph_peak_caller.peakscores import ScoredPeak
if pytest.__version__ < "3.0.0":
    pytest.skip()

@pytest.mark.skip("Legacy")
class TestMaxPath(unittest.TestCase):

    def test_find_max_path_through_subgraph_two_node_graph(self):

        graph = Graph({
                1: Block(10),
                2: Block(10)
            },
            {
                1: [2]
            }
        )

        peak = ConnectedAreas(graph,
                              {
                                  2: [0, 4],
                                  1: [5, 10]
                              })

        binary_peak = BinaryContinousAreas.from_old_areas(peak)
        qvalues = DensePileup.from_base_value(graph, 10)
        print("q values")
        print(qvalues)
        print(qvalues.data._values)
        scored_peak = ScoredPeak.from_peak_and_pileup(binary_peak, qvalues)
        print(scored_peak)

        max_path = scored_peak.get_max_path()

        self.assertEqual(max_path, Interval(5, 4, [1, 2]))

    def test_find_max_path_through_subgraph_multiple_paths(self):

        graph = Graph({
                1: Block(10),
                2: Block(10),
                3: Block(10),
                4: Block(10)
            },
            {
                1: [2, 3],
                2: [4],
                3: [4]
            }
        )

        peak = ConnectedAreas(graph,
                              {
                                  2: [0, 10],
                                  3: [0, 10],
                                  1: [5, 10],
                                  4: [0, 3]
                              })


        binary_peak = BinaryContinousAreas.from_old_areas(peak)
        qvalues = DensePileup.from_intervals(graph,
                [
                    Interval(7, 2, [1, 3, 4])  # Giving higher qvalue
                                                # through this path
                ])

        print(qvalues)

        scored_peak = ScoredPeak.from_peak_and_pileup(binary_peak, qvalues)
        print(scored_peak)

        max_path = scored_peak.get_max_path()
        self.assertEqual(max_path, Interval(5, 3, [1, 3, 4]))

    def test_find_max_path_on_start_and_end_node(self):

        graph = Graph({
                1: Block(10),
                2: Block(10),
                3: Block(10),
                4: Block(10)
            },
            {
                1: [2, 3],
                2: [4],
                3: [4]
            }
        )

        peak = ConnectedAreas(graph,
                              {
                                  2: [0, 10],
                                  4: [0, 10],
                              })

        binary_peak = BinaryContinousAreas.from_old_areas(peak)
        qvalues = DensePileup.from_intervals(
            graph,
            [
                Interval(7, 2, [1, 2, 4])
            ])
        scored_peak = ScoredPeak.from_peak_and_pileup(binary_peak, qvalues)

        max_path = scored_peak.get_max_path()
        self.assertEqual(max_path, Interval(0, 10, [2, 4]))

    def test_find_max_path_through_subgraph_with_illegal_paths(self):

        graph = Graph({
                1: Block(10),
                2: Block(10),
                3: Block(10),
                4: Block(10)
            },
            {
                1: [2, 3],
                2: [4],
                -4: [-3]   # Making 3=>4 not allowed path
            }
        )

        peak = ConnectedAreas(graph,
                              {
                                  2: [0, 10],
                                  3: [0, 10],
                                  1: [5, 10],
                                  4: [0, 8]
                              })

        binary_peak = BinaryContinousAreas.from_old_areas(peak)
        qvalues = DensePileup.from_intervals(graph,
                [
                    Interval(0, 10, [3]), # Higher value on 3 than 2
                    Interval(0, 10, [3]),
                    Interval(0, 10, [4]), # Highest value if ending on 4
                    Interval(0, 10, [4]),
                    Interval(0, 10, [1]), # Highest value if inncluding 1
                    Interval(0, 10, [1]), # Highest value if inncluding 1
                    Interval(0, 10, [1, 2, 4])
                ])

        scored_peak = ScoredPeak.from_peak_and_pileup(binary_peak, qvalues)

        max_path = scored_peak.get_max_path()
        print(max_path)

        self.assertEqual(max_path, Interval(5, 8, [1, 2, 4]))


if __name__ == "__main__":
    unittest.main()
