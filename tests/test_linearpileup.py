import unittest
from graph_peak_caller.linearsnarls import *


class TestLinearPileup(unittest.TestCase):

    def setUp(self):
        starts = np.array([1, 4, 6])
        ends = np.array([2, 5, 10])
        self.pileup = LinearPileup.create_from_starts_and_ends(starts, ends)

        self.valued_pileup = LinearPileup(np.array([1, 4, 8]), np.array([1, 0, 3]))

    def test_create_from_starts_and_ends(self):
        self.assertEqual(LinearPileup(np.array([1, 2, 4, 5, 6, 10]),
                                      np.array([1., 0., 1., 0., 1., 0.])),
                         self.pileup)

    def test_threshold(self):
        self.valued_pileup.threshold(2)
        self.assertTrue(np.all(self.valued_pileup.values == [2, 2, 3]))

    def test_maximum(self):
        pileup1 = LinearPileup(np.array([0,   3, 5.2, 10.8]),
                               np.array([10, 20,  15,   25]))

        pileup2 = LinearPileup(np.array([0,  5.2, 10]),
                               np.array([25,  10, 15]))
        pileup1.maximum(pileup2)

        true_max = LinearPileup(np.array([0, 5.2, 10.8]),
                                np.array([25, 15, 25]))
        self.assertEqual(pileup1, true_max)

    def _test_max(self):
        pileup1 = LinearPileup(np.array([1, 5, 10]), np.array([1, 2, 3]))
        pileup2 = LinearPileup(np.array([1, 5, 10]), np.array([3, 2, 1]))

        pileup1.maximum(pileup2)
        print(pileup1.values)
        self.assertTrue(np.all(pileup1.values == [3, 2, 3]))


class TestLinearIntervalCollection(object):

    def test_extend_non_overlapping(self):

        collection = LinearIntervalCollection(
            np.array([5, 10]),
            np.array([7, 12])
        )

        pileup = collection.extend(2)
        print(pileup.indices)
        self.assertTrue(np.all(pileup.indices == [4, 8, 9, 13]))
        self.assertTrue(np.all(pileup.values == [1, 0, 1, 0]))

    def test_extend_overlapping(self):

        collection = LinearIntervalCollection(
            np.array([5, 10]),
            np.array([7, 12])
        )

        pileup = collection.extend(5)
        print(pileup.indices)
        self.assertTrue(np.all(pileup.indices == [1, 6, 11, 16]))
        self.assertTrue(np.all(pileup.values == [1, 2, 1, 0]))

    def test_extend_overlapping_and_touching(self):

        collection = LinearIntervalCollection(
            np.array([0, 5, 8]),
            np.array([4, 7, 10])
        )

        pileup = collection.extend(2)
        print(pileup.indices)
        print(pileup.values)
        self.assertTrue(np.all(pileup.indices == [0, 7, 8, 11]))
        self.assertTrue(np.all(pileup.values == [1, 2, 1, 0]))


if __name__ == "__main__":
    unittest.main()
