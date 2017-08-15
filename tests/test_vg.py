import unittest
import json
from binary_correlation.vg import Position
import offsetbasedgraph
position_json = json.loads("{offset: 10, node_id: 0}")
position_json2 = json.loads("{offset: 0, node_id: 1, is_reverse: true}")


class TestPosition(unittest.TestCase):

    def test_json(self):
        position = Position.from_json(position_json)
        self.assertEqual(position.offset, 10)
        self.assertEqual(position.node_id, 0)
        self.assertEqual(position.is_reverse, False)

    def test_json(self):
        position = Position.from_json(position_json2)
        self.assertEqual(position.offset, 0)
        self.assertEqual(position.node_id, 1)
        self.assertEqual(position.is_reverse, True)

    def test_translate(self):
        position = Position(10, 20, False)
        obg_position = offsetbasedgraph.Position(10, 20)
        trans_position = position.to_obg()
        self.assertEqual(trans_position, obg_position)


class TestPath(unittest.TestCase):
    def test_path(self):
        pass
