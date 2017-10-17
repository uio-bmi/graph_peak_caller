import logging
logging.basicConfig(level=logging.ERROR)
from graph_peak_caller import callpeaks
from offsetbasedgraph import IntervalCollection
import offsetbasedgraph as obg
import cProfile
import pyvg
from pyvg.util import vg_gam_file_to_interval_collection
from pyvg.sequences import SequenceRetriever

import traceback
import warnings
import sys

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file, 'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = warn_with_traceback


def run_with_gam(gam_file_name, vg_graph_file_name,
                 limit_to_chromosomes=False):
    # logging.basicConfig(level=logging.error)
    vg_graph = pyvg.Graph.create_from_file(vg_graph_file_name)
    ob_graph = vg_graph.get_offset_based_graph()
    ob_graph.to_file("obgraph")
    # ob_graph = obg.GraphWithReversals.from_file("obgraph")

    print(ob_graph.adj_list[68566])
    print(ob_graph.adj_list[68567])

    # print(ob_graph.blocks)
    reads_intervals = vg_gam_file_to_interval_collection(
         None, gam_file_name, ob_graph, max_intervals=20000)

    # reads_intervals.to_file("test_obg_intervals")
    # exit()
    #IntervalCollection.interval_class = obg.DirectedInterval
    #reads_intervals = IntervalCollection.from_file("test_obg_intervals", graph=ob_graph)
    #control_intervals = IntervalCollection.from_file("test_obg_intervals", graph=ob_graph)
    control_intervals = vg_gam_file_to_interval_collection(
         None, gam_file_name, ob_graph, max_intervals=20000)

    experiment_info = callpeaks.ExperimentInfo(12000000, 103, 50)
    caller = callpeaks.CallPeaks(
        ob_graph, reads_intervals, control_intervals, experiment_info=experiment_info,
        out_file_base_name="real_data_", has_control=False)
    caller.verbose = True
    caller.run()
    retriever = SequenceRetriever.from_vg_graph("cactus-mhc.vg")
    sequences = [retriever.get_interval_sequence(max_path)
                 for max_path in caller.max_paths]
    f = open("real_data_sequences", "w")
    for seq in sequences:
        f.write(seq + "\n")

if __name__ == "__main__":
    #run_with_gam("reads3_large.gam", "../graph_peak_caller/dm_test_data/x.json", limit_to_chromosomes=["chr3R", "chr3L", "chr2R", "chrX", "chr2L", "chrY", "chr4"])
    dm_folder = "../graph_peak_caller/dm_test_data/"
    cProfile.run('run_with_gam("ENCFF291CUA.gam", "cactus-mhc.json")')
