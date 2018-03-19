import unittest
from offsetbasedgraph import GraphWithReversals, Block, \
    IntervalCollection, DirectedInterval as Interval, SequenceGraph
from graph_peak_caller.peakcollection import Peak, PeakCollection
import os
from graph_peak_caller.sparsediffs import SparseValues
from graph_peak_caller.command_line_interface import run_argument_parser
import numpy as np
from graph_peak_caller.peakfasta import PeakFasta


class Arguments(object):
    def __init__(self, dict):
        self.__dict__ = dict


class TestWrapper(unittest.TestCase):

    def setUp(self):
        self.correct_ob_graph = GraphWithReversals(
            {1: Block(7), 2: Block(4), 3: Block(7), 4: Block(4)},
            {1: [2, 3], 2: [4], 3: [4]})
        self.correct_ob_graph.convert_to_numpy_backend()

        self.correct_sequence_graph = SequenceGraph.create_empty_from_ob_graph(self.correct_ob_graph)
        self.correct_sequence_graph.set_sequences_using_vg_json_graph("tests/vg_test_graph.json")

        remove_files = ["tests/testgraph.obg", "tests/test_linear_map_starts.pickle",
                        "tests/test_linear_map_ends.pickle", "tests/test_linear_map.length",
                        "tests/sample.intervalcollection", "tests/testintervals.intervalcollection",
                        "tests/testsequences.fasta"]
        for file in remove_files:
            if os.path.isfile(file):
                os.remove(file)


class TestCommandLineInterface(TestWrapper):
    def test_create_ob_graph(self):
        run_argument_parser(["create_ob_graph", "-o", "tests/testgraph.obg",
                             "tests/vg_test_graph.json"])
        graph = GraphWithReversals.from_numpy_file("tests/testgraph.obg")
        self.assertEqual(graph, self.correct_ob_graph)

    def test_create_ob_graph_no_output_name(self):
        run_argument_parser(["create_ob_graph", "tests/vg_test_graph.json"])
        graph = GraphWithReversals.from_numpy_file("tests/vg_test_graph.nobg")
        self.assertEqual(graph, self.correct_ob_graph)

    def test_all_steps(self):
        run_argument_parser(["create_ob_graph", "-o",
                             "tests/testgraph.obg",
                             "tests/vg_test_graph.json"])
        run_argument_parser(['create_linear_map', "--graph",
                             "tests/testgraph.obg"])

        IntervalCollection([Interval(1, 1, [1, 2])]).to_file(
            "tests/sample.intervalcollection")

        run_argument_parser(["callpeaks",
                             "--graph", "tests/testgraph.obg",
                             "-s", "tests/sample.intervalcollection",
                             "-n", "tests/test_experiment_",
                             "-f", "10",
                             "-r", "7"])

    def _test_multigraph(self):
        run_argument_parser(["callpeaks_whole_genome",
                             "--chromosomes",
                             "--graph", "tests/testgraph.obg",
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

    def test_intervals_to_fasta_from_fasta(self):
        run_argument_parser(["create_ob_graph", "-o",
                             "tests/testgraph.obg",
                             "tests/vg_test_graph.json"])


        PeakCollection([Peak(0, 2, [1, 2], score=3)]).to_file("tests/testintervals.intervalcollection",
                                                     text_file=True)
        run_argument_parser(["peaks_to_fasta", "tests/testgraph.obg.sequences",
                            "tests/testintervals.intervalcollection", "tests/testsequences.fasta"])

        collection = PeakCollection.from_fasta_file("tests/testsequences.fasta")
        self.assertEqual(len(collection.intervals), 1)
        self.assertEqual(collection.intervals[0].sequence.lower(), "tttcccctt")

    def test_get_summits(self):

        qvalues = SparseValues(np.array([0]), np.array([3]))
        qvalues.track_size = 22
        qvalues.to_sparse_files("tests/test_qvalues")

        run_argument_parser(["create_ob_graph", "-o",
                             "tests/testgraph.obg",
                             "tests/vg_test_graph.json"])
        max_paths = PeakCollection([Peak(0, 2, [1, 2], score=3)])
        PeakFasta(self.correct_sequence_graph).write_max_path_sequences("tests/test_max_paths.fasta", max_paths)

        run_argument_parser(["get_summits", "-g", "tests/testgraph.obg",
                            "tests/test_max_paths.fasta", "tests/test_qvalues", "2"])

        result = PeakCollection.from_fasta_file("tests/test_max_paths_summits.fasta")
        self.assertEqual(result.intervals[0], Peak(2, 6, [1]))
        self.assertEqual(result.intervals[0].sequence.lower(), "tccc")


if __name__ == "__main__":
    unittest.main()



