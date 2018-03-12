import unittest
from graph_peak_caller.command_line_interface import create_ob_graph, \
    create_linear_map_interface, run_callpeaks_interface, concatenate_sequence_files
from offsetbasedgraph import GraphWithReversals, Block, IntervalCollection, DirectedInterval as Interval
from graph_peak_caller.peakcollection import Peak
import os
import sys
from unittest.mock import patch
from graph_peak_caller.command_line_interface import run_argument_parser

class Arguments(object):
    def __init__(self, dict):
        self.__dict__ = dict

class TestWrapper(unittest.TestCase):

    def setUp(self):
        self.correct_ob_graph = GraphWithReversals(
            {1: Block(7), 2: Block(4), 3: Block(7), 4: Block(4)},
            {1: [2, 3], 2: [4], 3: [4]})

        remove_files = ["tests/testgraph.obg", "tests/test_linear_map_starts.pickle",
                        "tests/test_linear_map_ends.pickle", "tests/test_linear_map.length",
                        "tests/sample.intervalcollection"]
        for file in remove_files:
            if os.path.isfile(file):
                os.remove(file)

class TestCommandLineInterface(TestWrapper):
    def test_create_ob_graph(self):
        run_argument_parser(["create_ob_graph" , "-o", "tests/testgraph.obg",
                             "tests/vg_test_graph.json"])
        graph = GraphWithReversals.from_numpy_file("tests/testgraph.obg")
        self.assertEqual(graph, self.correct_ob_graph)

    def test_create_ob_graph_no_output_name(self):
        run_argument_parser(["create_ob_graph", "tests/vg_test_graph.json"])
        graph = GraphWithReversals.from_numpy_file("tests/vg_test_graph.nobg")
        self.assertEqual(graph, self.correct_ob_graph)

    def test_all_steps(self):
        run_argument_parser(["create_ob_graph" , "-o", "tests/testgraph.obg",
                             "tests/vg_test_graph.json"])
        run_argument_parser(['create_linear_map', "--graph", "tests/testgraph.obg",
                             '-o', 'tests/test_linear_map',
                             'tests/vg_test_graph.snarls'])

        IntervalCollection([Interval(1, 1, [1, 2])]).to_file("tests/sample.intervalcollection")

        run_argument_parser(["callpeaks",
                             "--graph", "tests/testgraph.obg",
                             "-m", "tests/test_linear_map",
                             "-s", "tests/sample.intervalcollection",
                             "-n", "tests/test_experiment_",
                             "-f", "10",
                             "-r", "7"])

    def test_concatenate_fasta_files(self):
        if os.path.isfile("out.fasta"):
            os.remove("out.fasta")

        file1 = open("1_sequences.fasta", "w")
        peak1 = Peak(3, 5, [1], score=3)
        peak1_header = ">peak1 " + peak1.to_file_line() + "\n"
        file1.writelines([peak1_header])
        file1.writelines(["AACC\n"])
        file1.close()

        file2 = open("2_sequences.fasta", "w")
        peak2 = Peak(2, 5, [1], score=7)
        peak2_header = ">peak2 " + peak2.to_file_line() + "\n"
        file2.writelines([peak2_header])
        file2.writelines(["CCAA\n"])
        file2.close()

        run_argument_parser(["concatenate_sequence_files", "1,2", "out.fasta"])

        outfile = open("out.fasta")
        lines = list(outfile.readlines())
        self.assertEqual(Peak.from_file_line(lines[0].split(maxsplit=1)[1]), peak2)
        self.assertEqual(Peak.from_file_line(lines[2].split(maxsplit=1)[1]), peak1)
        outfile.close()


"""
class TestCommandLineInterface(TestWrapper):


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
                'graph': GraphWithReversals.from_file("tests/testgraph.obg"),
                "vg_snarls_file_name": "tests/vg_test_graph.snarls",
                "out_file_base_name": "tests/test_linear_map"
            }
        ))

        sample_reads = IntervalCollection([Interval(1, 1, [1, 2])])
        control_reads = IntervalCollection([Interval(1, 1, [1, 2])])

        run_callpeaks_interface(Arguments(
            {
                'graph': GraphWithReversals.from_file("tests/testgraph.obg"),
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

    def test_concatenate_fasta_files(self):
        if os.path.isfile("out.fasta"):
            os.remove("out.fasta")

        file1 = open("1_sequences.fasta", "w")
        peak1 = Peak(3, 5, [1], score=3)
        peak1_header = ">peak1 " + peak1.to_file_line() + "\n"
        file1.writelines([peak1_header])
        file1.writelines(["AACC\n"])
        file1.close()

        file2 = open("2_sequences.fasta", "w")
        peak2 = Peak(2, 5, [1], score=7)
        peak2_header = ">peak2 " + peak2.to_file_line() + "\n"
        file2.writelines([peak2_header])
        file2.writelines(["CCAA\n"])
        file2.close()

        args = Arguments({
            "chromosomes": "1,2",
            "out_file_name": "out.fasta"
        })

        concatenate_sequence_files(args)

        outfile = open("out.fasta")
        lines = list(outfile.readlines())
        self.assertEqual(Peak.from_file_line(lines[0].split(maxsplit=1)[1]), peak2)
        self.assertEqual(Peak.from_file_line(lines[2].split(maxsplit=1)[1]), peak1)
        outfile.close()
"""


if __name__ == "__main__":
    unittest.main()



