import unittest
from graph_peak_caller.sparsepileup import SparsePileup
from graph_peak_caller.pileupcleaner import PileupCleaner, IntervalWithinBlock
from cyclic_graph import get_small_cyclic_graph, get_large_cyclic_graph,\
    get_padded_cyclic_graph
from graph_peak_caller.pileupcleaner2 import PileupCleaner2
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
        filtered = cleaner.filter_on_length(5)
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
        filtered = cleaner.filter_on_length(16)
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
        filtered = cleaner.filter_on_length(4)
        print("Filtered")
        print(filtered)
        self.assertTrue(obg.Interval(0, 3, [1, 2]) in filtered)
        self.assertTrue(obg.Interval(0, 3, [2, 3]) in filtered)
        self.assertFalse(obg.Interval(0, 3, [1, 2, 3]) in filtered)

    def test_filter_on_length_long_interval(self):
        long_graph = obg.Graph(
            {
                1: obg.Block(10),
                2: obg.Block(10),
                3: obg.Block(10),
                4: obg.Block(10),
                5: obg.Block(10),
            },
            {1: [2],
             2: [3],
             3: [4],
             4: [5]}
        )
        pileup = SparsePileup.from_intervals(
            long_graph,
            [
                obg.Interval(0, 10, [2, 3, 4]),
            ]
        )
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        filtered = cleaner.filter_on_length(25)
        self.assertTrue(obg.Interval(0, 10, [2, 3, 4]) in filtered)

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
        filtered = cleaner.filter_on_length(11)
        print(filtered)
        self.assertTrue(obg.Interval(0, 3, [2, 5]) in filtered)
        self.assertFalse(obg.Interval(0, 3, [1, 5]) in filtered)


class TestCyclicCleanup(unittest.TestCase):
    def setUp(self):
        self.small_graph = get_small_cyclic_graph()
        self.large_graph = get_large_cyclic_graph()
        self.padded_graph = get_padded_cyclic_graph()

    def assertIntervalsGiveSamePileup(self, intervals, true_intervals):
        if not intervals:
            self.assertEqual(len(true_intervals), 0)

        pileup = SparsePileup.from_intervals(
            intervals[0].graph,
            intervals)
        pileup.threshold(1)
        true_pileup = SparsePileup.from_intervals(
            true_intervals[0].graph,
            true_intervals)
        true_pileup.threshold(1)
        self.assertEqual(pileup, true_pileup)

    def test_loop_with_surrounding(self):
        start_intervals = [obg.Interval(
            90, 100, [1, 2],
            graph=self.padded_graph)]
        pileup = SparsePileup.from_intervals(
            self.padded_graph, start_intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        intervals = cleaner.filter_on_length(400)
        self.assertIntervalsGiveSamePileup(intervals, start_intervals)

    def test_loop_with_surrounding_fail(self):
        start_intervals = [obg.Interval(
            90, 90, [1, 2],
            graph=self.padded_graph)]
        pileup = SparsePileup.from_intervals(
            self.padded_graph, start_intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        intervals = cleaner.filter_on_length(400)
        self.assertEqual(len(intervals), 0)

    def test_loop_with_start_and_end_intervals(self):
        start_intervals = [
            obg.Interval(
                90, 10, [1, 2],
                graph=self.padded_graph),
            obg.Interval(
                90, 100, [2],
                graph=self.padded_graph)]
        pileup = SparsePileup.from_intervals(
            self.padded_graph, start_intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner(pileup)
        intervals = cleaner.filter_on_length(400)
        self.assertEqual(intervals, [])


class TestPileupCleaner2(unittest.TestCase):

    def setUp(self):
        self.small_graph = get_small_cyclic_graph()
        self.large_graph = get_large_cyclic_graph()
        self.padded_graph = get_padded_cyclic_graph()

    def assert_interval_contained_in_intervals(self, interval, intervals):
        contained = False
        for other in intervals:
            if other.contains(interval):
                contained = True

        self.assertTrue(contained, "%s not in filtered" % str(interval))

    def assertIntervalsGiveSamePileup(self, intervals, true_intervals):
        if not intervals:
            self.assertEqual(len(true_intervals), 0)

        pileup = SparsePileup.from_intervals(
            intervals[0].graph,
            intervals)
        pileup.threshold(1)
        true_pileup = SparsePileup.from_intervals(
            true_intervals[0].graph,
            true_intervals)
        true_pileup.threshold(1)
        self.assertEqual(pileup, true_pileup)

    def test_create_interval_indices(self):
        graph = obg.Graph(
            {1: obg.Block(3),
             2: obg.Block(3),
             3: obg.Block(3)
             },
             {1: [2],
              -3: [-2]}
        )

        intervals = [
            obg.Interval(0, 3, [1]),
            obg.Interval(0, 3, [2]),
            obg.Interval(0, 3, [3])
            ]
        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        cleaner.create_interval_indices()

        start = cleaner.intervals_at_start_of_block
        print(start)
        end = cleaner.intervals_at_end_of_block
        self.assertTrue(intervals[0] in start[1])
        self.assertTrue(intervals[0] in end[1])

        self.assertTrue(intervals[1] in start[2])
        self.assertTrue(intervals[1] in end[2])

        self.assertTrue(intervals[2] in start[3])

    def test_filter_on_length_simple(self):
        graph = obg.Graph(
            {
            4: obg.Block(3),
            1: obg.Block(3),
             2: obg.Block(3),
             3: obg.Block(3)
             },
             {
                -2: [-4],
                1: [2],
                -3: [-2]}
        )
        intervals = [
            obg.Interval(2, 3, [4]),
            obg.Interval(0, 3, [1]),
            obg.Interval(0, 3, [2]),
            obg.Interval(1, 2, [3])
            ]
        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        filtered = cleaner.filter_on_length(1)

        self.assertTrue(obg.Interval(0, 3, [1, 2]) in filtered)
        self.assertTrue(obg.Interval(0, 1, [-2, -4]) in filtered)
        self.assertTrue(obg.Interval(1, 2, [3]) in filtered)

    def test_filter_on_length_simple2(self):
        graph = obg.Graph(
            {
                1: obg.Block(3),
                2: obg.Block(3),
                3: obg.Block(3),
                4: obg.Block(3)
            },
            {
                1: [2],
                2: [3],
                -3: [-2],
                -4: [-3]
            }
        )

        intervals = [
            obg.Interval(0, 3, [1, 2, 3]),
            obg.Interval(0, 3, [-4, -3, -2])
        ]
        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        filtered = cleaner.filter_on_length(7)
        self.assertTrue(obg.Interval(0, 3, [1, 2, 3]) in filtered)
        self.assertTrue(obg.Interval(0, 3, [-4, -3, -2]) in filtered)

    def test_filter_on_length_simple3(self):
        graph = obg.Graph(
            {
                1: obg.Block(3),
                2: obg.Block(3),
                3: obg.Block(3)
            },
            {
                1: [2],
                2: [3],
                -3: [-1]
            }
        )

        intervals = [
            obg.Interval(1, 3, [1, 2, 3])
        ]
        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        filtered = cleaner.filter_on_length(4)
        self.assertTrue(obg.Interval(1, 3, [1, 2, 3]) in filtered)
        self.assertTrue(obg.Interval(0, 2, [-3, -1]) in filtered)

    def test_filter_on_length_simple4(self):
        graph = obg.Graph(
            {
                1: obg.Block(3),
                2: obg.Block(3),
                3: obg.Block(3),
                4: obg.Block(3),
                5: obg.Block(3),
                6: obg.Block(3),
            },
            {
                1: [3],
                2: [3],
                3: [4],
                4: [5],
                5: [6]
            }
        )

        intervals = [
            obg.Interval(0, 3, [1]),
            obg.Interval(0, 3, [2]),
            obg.Interval(0, 3, [3, 4, 5, 6]),
        ]

        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        filtered = cleaner.filter_on_length(4)

        self.assertTrue(obg.Interval(0, 3, [1, 3, 4, 5, 6]) in filtered)
        self.assertTrue(obg.Interval(0, 3, [2, 3, 4, 5, 6]) in filtered)
        #self.assertTrue(obg.Interval(0, 2, [-3, -1]) in filtered)


    def test_filter_on_length_loop1(self):
        graph = obg.Graph(
            {
                1: obg.Block(3)
            },
            {
                1: [1]
            }
        )
        intervals = [
            obg.Interval(0, 3, [1])
        ]
        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        filtered = cleaner.filter_on_length(4)
        self.assertTrue(obg.Interval(0, 3, [1]) in filtered)


    def test_filter_on_length_loop2(self):
        graph = obg.Graph(
            {
                1: obg.Block(3),
                2: obg.Block(3)
            },
            {
                -2: [-1],
                -1: [-2]

            }
        )
        intervals = [
            obg.Interval(0, 3, [1, 2])
        ]
        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        filtered = cleaner.filter_on_length(4)
        self.assertTrue(obg.Interval(0, 3, [-2, -1]) in filtered)

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
        cleaner = PileupCleaner2(pileup)
        cleaner.find_trivial_intervals_within_blocks(cleaner.valued_areas)
        filtered = cleaner.filter_on_length(11)
        print(filtered)
        self.assertTrue(obg.Interval(0, 10, [1, 2, 3]) in filtered)
        self.assertTrue(obg.Interval(7, 10, [-5, -2]) in filtered)
        self.assertTrue(obg.Interval(0, 10, [-1, -3]) in filtered)
        self.assertFalse(obg.Interval(7, 10, [-5, -2, -1]) in filtered)

    def test_loop_with_surrounding(self):
        start_intervals = [obg.Interval(
            90, 100, [1, 2],
            graph=self.padded_graph)]
        pileup = SparsePileup.from_intervals(
            self.padded_graph, start_intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        intervals = cleaner.filter_on_length(400)
        self.assertIntervalsGiveSamePileup(intervals, start_intervals)

    def test_loop_with_surrounding_fail(self):
        start_intervals = [obg.Interval(
            90, 90, [1, 2],
            graph=self.padded_graph)]
        pileup = SparsePileup.from_intervals(
            self.padded_graph, start_intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        intervals = cleaner.filter_on_length(400)
        self.assertEqual(len(intervals), 0)

    def test_loop_with_start_and_end_intervals(self):
        start_intervals = [
            obg.Interval(
                90, 10, [1, 2],
                graph=self.padded_graph),
            obg.Interval(
                90, 100, [2],
                graph=self.padded_graph)]
        pileup = SparsePileup.from_intervals(
            self.padded_graph, start_intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        intervals = cleaner.filter_on_length(400)
        self.assertEqual(intervals, [])

    def test_loop_to_self(self):
        graph = obg.Graph(
            {1: obg.Block(3)},
            {1: [-1],
             -1: [1]}
        )
        intervals = [obg.Interval(0, 3, [1])]
        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        filtered = cleaner.filter_on_length(1)
        self.assertTrue(obg.Interval(0, 3, [1]) in filtered)
        self.assertTrue(obg.Interval(0, 3, [1, -1]) in filtered)
        self.assertFalse(obg.Interval(0, 3, [1, -1, 1]) in filtered)
        self.assertFalse(obg.Interval(0, 3, [1, -1, 1, -1]) in filtered)

    def test_loop_to_self_two_blocks(self):
        graph = obg.Graph(
            {1: obg.Block(3),
             2: obg.Block(3)},
            {1: [2],
             2: [-2],
             -2: [-1],
             -1: [1]}
        )
        intervals = [obg.Interval(0, 3, [1, 2])]
        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        filtered = cleaner.filter_on_length(1)
        #self.assertTrue(obg.Interval(0, 3, [1, 2]) in filtered)
        self.assertTrue(obg.Interval(0, 3, [1, 2, -2, -1]) in filtered)
        self.assertFalse(obg.Interval(0, 3, [1, 2, -2, -1, 1, 2]) in filtered)

    def test_double_loop_to_self_two_blocks(self):
        graph = obg.Graph(
            {1: obg.Block(3),
             2: obg.Block(3)},
            {1: [2],
             2: [-2, 2, -1],
             -2: [-1],
             -1: [1, 2],

             }
        )
        intervals = [obg.Interval(0, 3, [1, 2])]
        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        filtered = cleaner.filter_on_length(1)
        #self.assertTrue(obg.Interval(0, 3, [1, 2]) in filtered)
        #self.assert_interval_contained_in_intervals(obg.Interval(0, 3, [1, 2]), filtered)
        self.assertTrue(obg.Interval(0, 3, [1, 2, -2, -1]) in filtered)
        self.assertFalse(obg.Interval(0, 3, [1, 2, -2, -1, 1, 2]) in filtered)

    def test_loop_in_middle(self):
        graph = obg.Graph(
            {
                1: obg.Block(3),
                2: obg.Block(3),
                3: obg.Block(3),
             },
            {
                1: [2],
                2: [-2],
                2: [3]
            }
        )
        intervals = [obg.Interval(1, 2, [1, 2, 3])]
        pileup = SparsePileup.from_intervals(
            graph, intervals)
        pileup.threshold(1)
        cleaner = PileupCleaner2(pileup)
        filtered = cleaner.filter_on_length(1)

        self.assertTrue(obg.Interval(1, 2, [1, 2, 3]) in filtered)


from random import randrange, seed
class TestPileupCleaner2OnRandomGraphs(unittest.TestCase):

    def setUp(self):
        from collections import defaultdict

        n_blocks = 60 #  5
        n_edges = n_blocks + 60 # +2
        blocks = {}
        blocks_list = []
        for i in range(1, n_blocks + 1):
            blocks[i] = obg.Block(3)
            blocks_list.append(i)

        # Random edges
        edge_dict = defaultdict(list)
        for i in range(0, n_edges):
            start = blocks_list[randrange(0, len(blocks_list))]
            end = blocks_list[randrange(0, len(blocks_list))]

            if randrange(0, 2) == 1:
                start = -start

            if randrange(0, 2) == 1:
                end = -end

            if end == start or end == -start:
                continue

            if end not in edge_dict[start]:
                edge_dict[start].append(end)

        self.graph = obg.Graph(blocks, edge_dict)
        print(self.graph)
        # Create random interval
        start_rp = randrange(1, n_blocks + 1)
        start_offset = randrange(0, 3)
        end_offset = max(randrange(0, 3), start_offset + 1)
        rps = [start_rp]
        prev_rp = start_rp
        for i in range(0, n_blocks):
            edges = self.graph.adj_list[prev_rp]
            if len(edges) == 0:
                break
            rp = edges[randrange(0, len(edges))]
            prev_rp = rp
            if rp in rps:
                break

            rps.append(rp)

        self.interval = obg.Interval(start_offset, end_offset, rps, self.graph)

    def assertIntervalsGiveSamePileup(self, intervals, true_intervals):
        if not intervals:
            self.assertEqual(len(true_intervals), 0)
        print(" === Assert pileups equal ==")
        print(intervals)
        pileup = SparsePileup.from_intervals(
            intervals[0].graph,
            intervals)
        pileup.threshold(1)

        true_intervals_trivial = PileupCleaner2._intervals_to_trivial_intervals(true_intervals)
        true_pileup = SparsePileup.from_intervals(
            true_intervals[0].graph,
            true_intervals_trivial)
        true_pileup.threshold(1)
        print("Pileup")
        print(pileup)
        print("True pileup")
        print(true_pileup)
        self.assertEqual(pileup, true_pileup)

    def assert_interval_contained_in_intervals(self, interval, intervals):
        contained = False
        for other in intervals:
            if other.contains(interval):
                contained = True

        self.assertTrue(contained, "%s not in filtered" % str(interval))

    def assert_intervals_cover_same_areas(self, interval_set1, interval_set2):
        trivial1 = PileupCleaner2._intervals_to_trivial_intervals(interval_set1)
        trivial2 = PileupCleaner2._intervals_to_trivial_intervals(interval_set2)
        unique1 = PileupCleaner2._remove_interval_duplicates(trivial1)
        unique2 = PileupCleaner2._remove_interval_duplicates(trivial2)
        print("Unique1")
        print(unique1)
        print(unique2)
        self.assertEqual(len(unique1), len(unique2))

    def get_trivial_intervals(self):
        graph = obg.Graph(
            {1: obg.Block(3),
             2: obg.Block(3)},
            {1: [-2]}
        )
        interval = obg.Interval(2, 3, [1, -2], graph)
        trivial = PileupCleaner2._intervals_to_trivial_intervals([interval])
        unique = PileupCleaner2._remove_interval_duplicates(trivial)
        self.assertTrue(len(trivial) == 2)
        print("Trivial")
        print(trivial)
        self.assertTrue(len(unique) == 2)


    def test_random_graph(self):
        seed(1)
        for i in range(0, 200):
            print(" ======== TEst case =======")
            self.setUp()
            print(self.graph)
            print(self.interval)
            pileup = SparsePileup.from_intervals(
            self.graph, [self.interval])
            pileup.threshold(1)
            cleaner = PileupCleaner2(pileup)
            filtered = cleaner.filter_on_length(1, return_single_rp_intervals=True)

            #self.assert_interval_contained_in_intervals(self.interval, filtered)
            self.assert_intervals_cover_same_areas(filtered, [self.interval])

if __name__ == "__main__":
    unittest.main()
