import unittest
import json
from graph_peak_caller.vg import *
import offsetbasedgraph
position_jsons = [json.loads(pos_str) for pos_str in
                  ('{"offset": 10, "node_id": 0}',
                   '{"offset": 0, "node_id": 1, "is_reverse": true}')]
positions = [Position(0, 10, False), Position(1, 0, True)]

edit_jsons = [json.loads(edit_str) for edit_str in
              ['{"to_length": 4, "from_length": 4}',
               '{"to_length": 1, "from_length": 1, "sequence": "N"}']]
edits = [Edit(4, 4, None), Edit(1, 1, "N")]

mapping_jsons = [
    {"position": position_jsons[0], "edit": edit_jsons}]

mappings = [Mapping(positions[0], edits)]

path_jsons = [
    {"name": "path1",
     "mapping": mapping_jsons}]
paths = [Path("path1", mappings)]
intervals = [offsetbasedgraph.Interval(10, 15, [0])]


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


class TestPath(unittest.TestCase):
    def test_json(self):
        path = Path.from_json(path_jsons[0])
        self.assertEquals(path, paths[0])

    def test_obg(self):
        for path, interval in zip(paths, intervals):
            self.assertEqual(path.to_obg(), interval)


if __name__ == "__main__":
    unittest.main()
