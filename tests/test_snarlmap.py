import numpy as np
import unittest
import offsetbasedgraph as obg
from test_snarls import snarl_graph2
from graph_peak_caller.linearsnarls import UnmappedIndices, LinearPileup,\
    create_control_from_objs
from graph_peak_caller.snarlmaps import LinearSnarlMap
from graph_peak_caller.sparsepileup import ValuedIndexes, SparsePileup

graph = obg.GraphWithReversals(
    {3: obg.Block(20), 5: obg.Block(10),
     12: obg.Block(20), 13: obg.Block(21),
     }, {})


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

    def test_create_control(self):
        intervals = [obg.DirectedInterval(0, 20, [3]),
                     obg.DirectedInterval(0, 10, [5]),
                     obg.DirectedInterval(0, 21, [13])]
        mapped_intervals = self.snarl_map.map_interval_collection(
            intervals)
        linear_pileup = LinearPileup.create_from_starts_and_ends(
            mapped_intervals.starts,
            mapped_intervals.ends)
        valued_indexes = linear_pileup.to_valued_indexes(self.snarl_map)
        graph_pileup = SparsePileup(graph)
        graph_pileup.data = valued_indexes
        true_sparse_pileup = SparsePileup(graph)
        true_sparse_pileup.data = {3: ValuedIndexes([], [], 2, 20),
                                   12: ValuedIndexes([], [], 2, 20),
                                   13: ValuedIndexes([], [], 2, 21),
                                   5: ValuedIndexes([], [], 2, 10)}
        print(graph_pileup)
        self.assertEqual(graph_pileup, true_sparse_pileup)

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
        mapped_vis = self.snarl_map.to_graph_pileup(unmapped_indices)
        self.assertEqual(mapped_vis, vis)


class TestLinearPileupMap(TestSnarlMap):
    def test_to_valued_indexes(self):
        """[0, 5, 10, 15, 20, 25, 30]"""
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


if __name__ == "__main__":
    unittest.main()
