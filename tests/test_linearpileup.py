import unittest
from graph_peak_caller.control.linearsnarls import *
from graph_peak_caller.control.linearintervals import LinearIntervalCollection


class TestLinearPileup(unittest.TestCase):

    def setUp(self):
        starts = np.array([1, 4, 6])
        ends = np.array([2, 5, 10])
        self.pileup = LinearPileup.create_from_starts_and_ends(starts, ends)

        self.valued_pileup = LinearPileup(np.array([1, 4, 8]),
                                          np.array([1, 0, 3]))

    def test_create_from_starts_and_ends(self):
        self.assertEqual(LinearPileup(np.array([1, 2, 4, 5, 6, 10]),
                                      np.array([1., 0., 1., 0., 1., 0.])),
                         self.pileup)

    def test_maximum(self):
        pileup1 = LinearPileup(np.array([0,   3, 5.2, 10.8]),
                               np.array([10, 20,  15,   25]))

        pileup2 = LinearPileup(np.array([0,  5.2, 10]),
                               np.array([25,  10, 15]))
        pileup1.maximum(pileup2)

        true_max = LinearPileup(np.array([0, 5.2, 10.8]),
                                np.array([25, 15, 25]))
        self.assertEqual(pileup1, true_max)


class TestExtendLinearIntervalCollection(unittest.TestCase):

    def test_extend_non_overlapping(self):

        collection = LinearIntervalCollection(
            np.array([5, 10]),
            np.array([7, 12])
        )

        extended = collection.extend(2)
        pileup = LinearPileup.create_from_starts_and_ends(
            extended.starts,
            extended.ends
        )
        self.assertTrue(np.all(pileup.indices == [3, 7, 8, 12]))
        self.assertTrue(np.all(pileup.values == [1, 0, 1, 0]))

    def test_extend_overlapping(self):

        collection = LinearIntervalCollection(
            np.array([5, 10]),
            np.array([7, 12])
        )

        extended = collection.extend(5)
        pileup = LinearPileup.create_from_starts_and_ends(
            extended.starts,
            extended.ends
        )
        print(pileup.indices)
        self.assertTrue(np.all(pileup.indices == [0, 5, 10, 15]))
        self.assertTrue(np.all(pileup.values == [1, 2, 1, 0]))

    def test_extend_complex(self):

        collection = LinearIntervalCollection(
            np.array([10, 15, 18]),
            np.array([14, 17, 20])
        )

        extended = collection.extend(6)
        pileup = LinearPileup.create_from_starts_and_ends(
            extended.starts,
            extended.ends
        )

        self.assertTrue(np.all(pileup.indices == [4, 9, 12, 16, 21, 24]))
        self.assertTrue(np.all(pileup.values == [1, 2, 3, 2, 1, 0]))


if __name__ == "__main__":
    unittest.main()
