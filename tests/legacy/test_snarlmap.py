import pytest
if pytest.__version__ < "3.0.0":
    pytest.skip()

import numpy as np
import unittest
import offsetbasedgraph as obg
# from test_snarls import snarl_graph2
# from graph_peak_caller.control.linearsnarls import \
#     UnmappedIndices, LinearPileup
# from graph_peak_caller.control.snarlmaps import LinearSnarlMap

graph = obg.GraphWithReversals(
    {3: obg.Block(20), 5: obg.Block(10),
     12: obg.Block(20), 13: obg.Block(21),
     }, {})


@pytest.mark.skip()
class TestSnarlMap(unittest.TestCase):
    def setUp(self):
        self.snarl_map = LinearSnarlMap.from_snarl_graph(snarl_graph2, graph)
        self.graph_positions = [obg.Position(5, 4),
                                obg.Position(3, 4),
                                obg.Position(12, 4),
                                obg.Position(13, 4)]

        self.linear_positions = [4, 31/20*4, 10+21/20*4, 14]
        self.linear_positions = [p for p in self.linear_positions]
        self.graph_interval = obg.DirectedInterval(self.graph_positions[0],
                                                   self.graph_positions[2])

    # TODO: Is test wrong?
    def _test_create_control(self):
        intervals = [obg.DirectedInterval(0, 20, [3]),
                     obg.DirectedInterval(0, 10, [5]),
                     obg.DirectedInterval(0, 21, [13])]
        mapped_intervals = self.snarl_map.map_interval_collection(
            intervals)
        linear_pileup = LinearPileup.create_from_starts_and_ends(
            mapped_intervals.starts,
            mapped_intervals.ends)
        graph_pileup = linear_pileup.to_sparse_pileup(self.snarl_map)

        true_sparse_pileup = OldSparsePileup(graph)
        true_data = {3: ValuedIndexes([], [], 2, 20),
                       12: ValuedIndexes([], [], 2, 20),
                       13: ValuedIndexes([], [], 2, 21),
                       5: ValuedIndexes([], [], 2, 10)}
        true_sparse_pileup.data = OldSparsePileupData([(key, val) for key, val in true_data.items()], graph=graph)
        #print(true_sparse_pileup.data)

        print("Graph pileup")
        print(graph_pileup)
        print("True")
        print(true_sparse_pileup)

        self.assertTrue(graph_pileup.equals_old_sparse_pileup(true_sparse_pileup))

    def test_graph_position_to_linear(self):
        for graph_pos, lin_pos in zip(self.graph_positions,
                                      self.linear_positions):
            mapped_pos = self.snarl_map.graph_position_to_linear(graph_pos)
            self.assertEqual(mapped_pos, lin_pos)

    def test_map_graph_interval(self):
        mapped_interval = self.snarl_map.map_graph_interval(
            self.graph_interval)
        self.assertEqual(mapped_interval,
                         (self.linear_positions[0],
                          self.linear_positions[2]), [5, 12])

    def test_to_graph_pileup(self):
        """[0, 5, 10, 15, 20, 25, 30]"""
        all_indices = [0, 5, 10, 15, 20, 25, 30]
        values = {idx: i for i, idx in enumerate(all_indices)}
        unmapped_indices = {3: all_indices[:],
                            5: [0, 5],
                            12: all_indices[2:],
                            13: all_indices[2:]}

        unmapped_indices = {node_id: UnmappedIndices(
            idxs, [values[idx] for idx in idxs])
                            for node_id, idxs in unmapped_indices.items()}

        vis = {3: (np.array(all_indices)*20/31, np.array(list(range(7)))),
               5: (np.array([0, 5]), np.array([0, 1])),
               12: ((np.array(all_indices[2:])-10)*20/21,
                    np.array(list(range(2, 7)))),
               13: ((np.array(all_indices[2:])-10,
                     np.array(list(range(2, 7)))))}
        vis = {node_id: ValuedIndexes(
            val[0][1:].astype("int"), val[1][1:],
            val[1][0], graph.node_size(node_id))
               for node_id, val in vis.items()}

        correct_pileup = OldSparsePileup(graph)

        print("Old sparse pileup")
        print(correct_pileup)

        correct_pileup.data = OldSparsePileupData([(key, val) for key, val in vis.items()], graph=graph)
        print("Correct pileup data")
        print(correct_pileup.data)
        print("Correct pileup")
        print(correct_pileup)
        correct_pileup = DensePileup.create_from_old_sparsepileup(correct_pileup)

        pileup = self.snarl_map.to_dense_pileup(unmapped_indices)

        print("Pileup from test")
        print(pileup)

        print("Correct pileup")
        print(correct_pileup)

        self.assertEqual(pileup, correct_pileup)


# Needs to be rewritten. To valued indexes does not exist anymore
"""
class TestLinearPileupMap(TestSnarlMap):
    def _test_to_valued_indexes(self):
        #[0, 5, 10, 15, 20, 25, 30]
        all_indices = [0, 5, 10, 15, 20, 25, 30]
        linear_pileup = LinearPileup(np.array(all_indices),
                                     np.array(list(range(7))))

        vis = {3: (np.array(all_indices)*20/31, np.array(list(range(7)))),
               5: (np.array([0, 5]), np.array([0, 1])),
               12: ((np.array(all_indices[2:])-10)*20/21,
                    np.array(list(range(2, 7)))),
               13: ((np.array(all_indices[2:])-10,
                     np.array(list(range(2, 7)))))}
        vis = {node_id: ValuedIndexes(
            val[0][1:].astype("int"), val[1][1:], val[1][0],
            graph.node_size(node_id))
               for node_id, val in vis.items()}
        mapped_vis = linear_pileup.to_valued_indexes(self.snarl_map)
        self.assertEqual(mapped_vis, vis)
"""

if __name__ == "__main__":
    unittest.main()
