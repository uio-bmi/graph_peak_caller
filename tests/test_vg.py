import unittest
from examples import *
import offsetbasedgraph


class TestPosition(unittest.TestCase):

    def test_json(self):
        for position_json, true_position in zip(position_jsons, positions):
            position = Position.from_json(position_json)
            self.assertEqual(position, true_position)

    def test_translate(self):
        position = Position(10, 20, False)
        obg_position = offsetbasedgraph.Position(10, 20)
        trans_position = position.to_obg()
        self.assertEqual(trans_position, obg_position)


class TestMapping(unittest.TestCase):
    def test_json(self):
        for mapping_json, true_mapping in zip(mapping_jsons, mappings):
            mapping = Mapping.from_json(mapping_json)
            self.assertEqual(mapping, true_mapping)

    def test_is_reverse(self):
        self.assertFalse(mappings[0].is_reverse())
        self.assertTrue(mappings[1].is_reverse())


class TestPath(unittest.TestCase):
    def test_json(self):
        path = Path.from_json(path_jsons[0])
        self.assertEquals(path, paths[0])

    def _est_obg(self):
        for path, interval in zip(paths, intervals):
            self.assertEqual(path.to_obg(), interval)


if __name__ == "__main__":
    unittest.main()
