import unittest
from graph_peak_caller.command_line_interface import create_ob_graph
from offsetbasedgraph import GraphWithReversals, Block

class Arguments(object):
    def __init__(self, dict):
        self.__dict__ = dict


vg_json_graph = \
"""{
"node": [
    {"id": 1, "sequence": "TTTCCCC"},
    {"id": 2, "sequence": "TTTT"},
    {"id": 3, "sequence": "CCCCTTT"}
],
"edge": [
    {"from": 1, "to": 2},
    {"from": 2, "to": 3}
]
}"""

class TestCommandLineInterface(unittest.TestCase):

    def setUp(self):
        with open("vg_graph.json", "w") as f:
            f.write(vg_json_graph.replace("\n", " "))

    def test_create_ob_graph(self):
        create_ob_graph(Arguments(
            {
                "vg_json_file_name": "vg_graph.json",
                "out_file_name": "testgraph.obg"
            })
        )
        graph = GraphWithReversals.from_numpy_file("testgraph.obg")
        self.assertEqual(graph, GraphWithReversals(
            {1: Block(7), 2: Block(4), 3: Block(7)},
            {1: [2], 2: [3]}))




if __name__ == "__main__":
    unittest.main()



