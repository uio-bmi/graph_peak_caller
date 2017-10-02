import unittest
from graph_peak_caller.sparsepileup import SparsePileup
from graph_peak_caller.pileupcleaner import PileupCleaner,  IntervalWithinBlock
import offsetbasedgraph as obg

class TestPileupCleaner(unittest.TestCase):

    def setUp(self):
        self.graph = obg.Graph(
            {1: obg.Block(10),
             2: obg.Block(10),
             3: obg.Block(10),
             4: obg.Block(10)}
            ,
            {1: [2, 4],
             2: [3],
             4: [3]
             }
        )
        #self.trivial_pileup = SparsePileup(self.graph)

        self.trivial_pileup = SparsePileup.from_intervals(
            self.graph,
            [
            obg.Interval(0, 5, [1]),
            obg.Interval(3, 5, [2]),
            obg.Interval(3, 10, [3]),
            obg.Interval(4, 6, [1])
            ]
        )
        self.trivial_pileup.threshold(1)

    def test_find_maximally_expanded_holes_one_block(self):
        graph = obg.Graph({
            1: obg.Block(10)
        }, {})
        pileup =  SparsePileup.from_intervals(
            graph,
            [
                obg.Interval(4, 5, [1], graph)
            ]
        )
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        holes = cleaner.find_maximally_expanded_holes()
        self.assertTrue(obg.Interval(0, 4, [1]) in holes)
        self.assertTrue(obg.Interval(5, 10, [1]) in holes)
        self.assertEqual(len(holes), 2)


    def test_find_maximally_expanded_holes_two_blocks(self):
        graph = obg.Graph({
                1: obg.Block(10),
                2: obg.Block(10)
            },
            {
                1: [2]
            }
        )
        pileup =  SparsePileup.from_intervals(
            graph,
            [
                obg.Interval(4, 5, [1], graph)
            ]
        )
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        holes = cleaner.find_maximally_expanded_holes()
        self.assertTrue(obg.Interval(0, 4, [1]) in holes)
        self.assertTrue(obg.Interval(5, 10, [1, 2]) in holes)
        self.assertFalse(obg.Interval(5, 10, [1]) in holes)
        self.assertFalse(obg.Interval(1, 10, [2]) in holes)


    def test_find_maximally_expanded_holes_complex(self):
        graph = obg.Graph({
                1: obg.Block(10),
                2: obg.Block(10),
                3: obg.Block(10),
                4: obg.Block(10),
                5: obg.Block(10)
            },
            {
                1: [3],
                2: [3],
                3: [4, 5]
            }
        )
        pileup =  SparsePileup.from_intervals(
            graph,
            [
                obg.Interval(0, 5, [1], graph),
                obg.Interval(5, 10, [2], graph),
                obg.Interval(3, 6, [3], graph),
                obg.Interval(5, 10, [4], graph),
                obg.Interval(0, 5, [5], graph),
            ]
        )
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        holes = cleaner.find_maximally_expanded_holes()
        print(holes)
        self.assertTrue(obg.Interval(0, 5, [2]) in holes)
        self.assertTrue(obg.Interval(0, 3, [3]) in holes)
        self.assertTrue(obg.Interval(6, 10, [3]) in holes)
        self.assertTrue(obg.Interval(6, 5, [3, 4]) in holes)
        self.assertFalse(obg.Interval(0, 5, [4]) in holes)
        self.assertFalse(obg.Interval(5, 10, [1]) in holes)

    def test_find_maximall_expanded_holes_cyclic(self):
        graph = obg.Graph({
                1: obg.Block(10),
                2: obg.Block(10),
                3: obg.Block(10),
                4: obg.Block(10),
                5: obg.Block(10)
            },
            {
                1: [3],
                2: [3],
                3: [4, 5, 1],
                -4: [2]
            }
        )
        pileup =  SparsePileup.from_intervals(
            graph,
            [
                obg.Interval(1, 5, [1], graph),
                obg.Interval(0, 10, [2], graph),
                obg.Interval(3, 6, [3], graph),
                obg.Interval(5, 10, [4], graph),
                obg.Interval(0, 5, [5], graph),
            ]
        )
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        holes = cleaner.find_maximally_expanded_holes()
        print(holes)
        self.assertFalse(obg.Interval(0, 5, [2]) in holes)
        self.assertTrue(obg.Interval(0, 3, [3]) in holes)
        self.assertTrue(obg.Interval(6, 10, [3]) in holes)
        self.assertTrue(obg.Interval(6, 5, [3, 4]) in holes)
        self.assertTrue(obg.Interval(6, 1, [3, 1]) in holes)
        self.assertTrue(obg.Interval(0, 5, [4]) in holes)
        self.assertFalse(obg.Interval(5, 10, [1]) in holes)

    def test_find_trivial_intervals(self):

        cleaner = PileupCleaner(self.trivial_pileup)
        trivial_intervals = cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)

        self.assertEqual(len(trivial_intervals), 3)

        self.assertTrue(obg.Interval(0, 6, [1]) in trivial_intervals)
        self.assertTrue(obg.Interval(3, 5, [2]) in trivial_intervals)
        self.assertTrue(obg.Interval(3, 10, [3]) in trivial_intervals)

    def test_IntervalWithinBlock(self):

        interval = IntervalWithinBlock(-1, 0, 5, [1], self.graph)
        self.assertTrue(interval.is_at_beginning_of_block())
        self.assertFalse(interval.is_at_end_of_block())
        self.assertTrue(len(interval.blocks_going_into()) == 0)
        self.assertTrue(len(interval.blocks_going_out_from()) == 0)

        interval = IntervalWithinBlock(-1, 0, 10, [2], self.graph)
        print(" === CASE === ")
        print(self.graph.reverse_adj_list)
        self.assertTrue(interval.is_at_beginning_of_block())
        self.assertTrue(interval.is_at_end_of_block())
        self.assertEqual(interval.blocks_going_into(), [1])
        self.assertEqual(interval.blocks_going_out_from(), [3])

    def test_blocks_going_out_from_advanced(self):
        graph = obg.Graph(
            {1: obg.Block(1),
             2: obg.Block(1),
             3: obg.Block(1)
            },
            {
                1: [2],
                -3: [-1]
            }
        )
        print("== test_blocks_going_out_from_advanced ==")
        print(graph.adj_list)
        print(graph.reverse_adj_list)
        interval = IntervalWithinBlock(1, 0, 1, [1], graph)
        self.assertTrue(interval.blocks_going_out_from(), [1, 3])


    def test_merge_right(self):
        interval = IntervalWithinBlock(1, 5, 10, [1], self.graph)
        interval2 = IntervalWithinBlock(2, 0, 5, [2], self.graph)

        merged = interval.merge_right(interval2)
        self.assertEqual(merged, obg.Interval(5, 5, [1, 2], self.graph))

    def test_create_interval_indices(self):
        cleaner = PileupCleaner(self.trivial_pileup)

        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        cleaner.create_interval_indices()


        index = cleaner.intervals_at_start_of_block
        print(index)
        self.assertTrue(1 in index)
        self.assertTrue(index[1] == [obg.Interval(0, 6, [1])])

    def test_filter_on_length_trivial(self):
        cleaner = PileupCleaner(self.trivial_pileup)
        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        filtered = cleaner.filter_on_length(2)

        self.assertTrue(obg.Interval(0, 6, [1]) in filtered)
        self.assertTrue(obg.Interval(3, 5, [2]) in filtered)
        self.assertTrue(obg.Interval(3, 10, [3]) in filtered)

        filtered = cleaner.filter_on_length(3)
        self.assertTrue(obg.Interval(0, 6, [1]) in filtered)
        self.assertTrue(obg.Interval(3, 10, [3]) in filtered)


        filtered = cleaner.filter_on_length(6)
        self.assertTrue(obg.Interval(3, 10, [3]) in filtered)

    def test_filter_on_length_interval_touching_end(self):
        self.trivial_pileup = SparsePileup.from_intervals(
            self.graph,
            [
            obg.Interval(0, 10, [1]),
            obg.Interval(0, 3, [2]),
            obg.Interval(3, 10, [3]),
            obg.Interval(0, 3, [4])
            ]
        )
        self.trivial_pileup.threshold(1)
        cleaner = PileupCleaner(self.trivial_pileup)
        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        filtered = cleaner.filter_on_length(2)
        self.assertTrue(obg.Interval(3, 10, [3]) in filtered)
        self.assertTrue(obg.Interval(0, 3, [1, 2]) in filtered)
        self.assertTrue(obg.Interval(0, 3, [1, 4]) in filtered)
        #self.assertEqual(len(filtered), 3)

    def test_filter_on_length_interval_spanning_multiple_rps(self):
        self.trivial_pileup = SparsePileup.from_intervals(
            self.graph,
            [
            obg.Interval(4, 10, [1]),
            obg.Interval(0, 10, [2]),
            obg.Interval(0, 5, [3]),
            obg.Interval(0, 10, [4])
            ]
        )
        self.trivial_pileup.threshold(1)
        cleaner = PileupCleaner(self.trivial_pileup)
        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        filtered = cleaner.filter_on_length(2)
        print(filtered)
        self.assertTrue(obg.Interval(4, 5, [1, 2, 3]) in filtered)
        self.assertTrue(obg.Interval(4, 5, [1, 4, 3]) in filtered)
        #self.assertEqual(len(filtered), 2)

    def test_filter_on_length_simple_loop(self):
        loop_graph = obg.Graph({1: obg.Block(10)},
                               {1: [1]})
        pileup = SparsePileup.from_intervals(
            loop_graph,
            [
                obg.Interval(0, 10, [1])
            ]
        )
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        filtered = cleaner.filter_on_length(2)
        print("Filtered loop")
        print(filtered)
        self.assertEqual(filtered, [obg.Interval(0, 10, [1])])


    def test_filter_on_length_simple_loop(self):
        loop_graph = obg.Graph({1: obg.Block(10)},
                               {1: [1]})
        pileup = SparsePileup.from_intervals(
            loop_graph,
            [
                obg.Interval(0, 10, [1])
            ]
        )
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        filtered = cleaner.filter_on_length(2)
        print("Filtered loop")
        print(filtered)
        self.assertEqual(filtered, [obg.Interval(0, 10, [1])])

    def test_filter_multi_directed(self):
        graph = obg.Graph(
            {1: obg.Block(3),
             2: obg.Block(3),
             3: obg.Block(3),
             },
            {2: [3],
            -2: [-1]}
        )
        pileup = SparsePileup.from_intervals(
            graph,
            [
                obg.Interval(0, 3, [1]),
                obg.Interval(0, 3, [2]),
                obg.Interval(0, 3, [3])
            ]
        )
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        filtered = cleaner.filter_on_length(1)
        print("Filtered")
        print(filtered)
        self.assertTrue(obg.Interval(0, 3, [1, 2]) in filtered)
        self.assertTrue(obg.Interval(0, 3, [2, 3]) in filtered)
        self.assertFalse(obg.Interval(0, 3, [1, 2, 3]) in filtered)

    def test_filter_advanced_multi_directed(self):
        graph = obg.Graph(
            {
                1: obg.Block(10),
                2: obg.Block(10),
                3: obg.Block(10),
                4: obg.Block(10),
                5: obg.Block(10)
            },
            {
                1: [2, 4],
                2: [3],
                -1: [-3],
                -5: [-2]
            }
        )
        pileup = SparsePileup.from_intervals(
            graph,
            [
                obg.Interval(0, 10, [1, 2, 3]),
                obg.Interval(0, 3, [5])
            ]
        )
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        filtered = cleaner.filter_on_length(1)
        print(filtered)
        self.assertTrue(obg.Interval(0, 3, [2, 5]) in filtered)
        self.assertFalse(obg.Interval(0, 3, [1, 5]) in filtered)

    def test_return_pileup(self):
        cleaner = PileupCleaner(self.trivial_pileup)
        cleaner.filter_on_length_and_return_pileup(3)


if __name__ == "__main__":
    unittest.main()
