import unittest
import numpy as np
import offsetbasedgraph as obg
from graph_peak_caller.areas import ValuedAreas,\
    BinaryContinousAreas


class TestValuedAreas(unittest.TestCase):
    def test_add_binary_areas(self):
        touched_nodes = set()
        graph = obg.Graph({1: obg.Block(100)}, {})
        starts = np.array([0] + [i*10 for i in range(9)])
        ends = np.array([100] + [(i+2)*10 for i in range(9)])
        bin_areas = []
        for start, end in zip(starts, ends):
            binary_connected_areas = BinaryContinousAreas(graph)
            binary_connected_areas.add(1, start, end)
            binary_connected_areas.sanitize()
            bin_areas.append(binary_connected_areas)
        valued_areas = ValuedAreas(graph)
        for bin_area in bin_areas:
            valued_areas.add_binary_areas(bin_area, touched_nodes=touched_nodes)

        self.assertTrue(np.all(
            valued_areas.get_starts_array(1) == np.sort(starts)))

        self.assertTrue(np.all(
            valued_areas.get_ends_array(1) == np.sort(ends)))

    def test_add_pos_and_negbinary_areas(self):
        touched_nodes = set()
        graph = obg.Graph({1: obg.Block(100)}, {})
        pos_starts = np.array([0] + [i*10 for i in range(9)])
        pos_ends = np.array([100] + [(i+2)*10 for i in range(9)])
        neg_starts = np.array([0] + [i*10+5 for i in range(9)])
        neg_ends = np.array([100] + [(i+1)*10+5 for i in range(9)])
        bin_areas = []
        for start, end in zip(pos_starts, pos_ends):
            binary_connected_areas = BinaryContinousAreas(graph)
            binary_connected_areas.add(1, start, end)
            binary_connected_areas.sanitize()
            bin_areas.append(binary_connected_areas)
        for start, end in zip(neg_starts, neg_ends):
            binary_connected_areas = BinaryContinousAreas(graph)
            binary_connected_areas.add(-1, start, end)
            binary_connected_areas.sanitize()
            bin_areas.append(binary_connected_areas)
        valued_areas = ValuedAreas(graph)
        for bin_area in bin_areas:
            valued_areas.add_binary_areas(bin_area, touched_nodes=touched_nodes)

        true_starts = np.sort(np.concatenate(
            [pos_starts, 100-neg_ends]))
        true_ends = np.sort(np.concatenate(
            [pos_ends, 100-neg_starts]))
        self.assertTrue(np.all(
            np.sort(valued_areas.get_starts_array(1)) == true_starts))
        self.assertTrue(np.all(
            np.sort(valued_areas.get_ends_array(1)) == true_ends))


if __name__ == "__main__":
    unittest.main()
