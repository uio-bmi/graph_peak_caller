import unittest
from graph_peak_caller.command_line_interface import create_ob_graph, \
    create_linear_map_interface, run_callpeaks_interface
from offsetbasedgraph import GraphWithReversals, Block, IntervalCollection, DirectedInterval as Interval
import os

class Arguments(object):
    def __init__(self, dict):
        self.__dict__ = dict


class TestCommandLineInterface(unittest.TestCase):

    def setUp(self):
        remove_files = ["tests/testgraph.obg", "tests/test_linear_map_starts.pickle",
                        "tests/test_linear_map_ends.pickle", "tests/test_linear_map.length"]
        for file in remove_files:
            if os.path.isfile(file):
                os.remove(file)

    def test_create_ob_graph(self):
        create_ob_graph(Arguments(
            {
                "vg_json_file_name": "tests/vg_test_graph.json",
                "out_file_name": "tests/testgraph.obg"
            })
        )
        graph = GraphWithReversals.from_numpy_file("tests/testgraph.obg")
        self.assertEqual(graph, GraphWithReversals(
            {1: Block(7), 2: Block(4), 3: Block(7), 4: Block(4)},
            {1: [2, 3], 2: [4], 3: [4]}))


    def test_all_steps(self):
        create_ob_graph(Arguments(
            {
                "vg_json_file_name": "tests/vg_test_graph.json",
                "out_file_name": "tests/testgraph.obg"
            })
        )
        create_linear_map_interface(Arguments(
            {
                "obg_file_name": "tests/testgraph.obg",
                "vg_snarls_file_name": "tests/vg_test_graph.snarls",
                "out_file_base_name": "tests/test_linear_map"
            }
        ))

        sample_reads = IntervalCollection([Interval(1, 1, [1, 2])])
        control_reads = IntervalCollection([Interval(1, 1, [1, 2])])

        run_callpeaks_interface(Arguments(
            {
                'graph_file_name': "tests/testgraph.obg",
                'vg_graph_file_name': "tests/vg_test_graph.vg",
                'linear_map_base_name': "tests/test_linear_map",
                'sample_reads_file_name': sample_reads,
                'control_reads_file_name': control_reads,
                'with_control': False,
                'out_base_name': 'tests/test_experiment_',
                'fragment_length': 10,
                'read_length': 7
            }
        ))



if __name__ == "__main__":
    unittest.main()



