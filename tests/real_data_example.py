import logging
logging.basicConfig(level=logging.INFO)
from graph_peak_caller import callpeaks
from offsetbasedgraph import IntervalCollection
import offsetbasedgraph as obg
import cProfile
import pyvg
from pyvg.util import vg_gam_file_to_interval_collection
from pyvg.sequences import SequenceRetriever
from graph_peak_caller.sparsepileup import SparsePileup
from graph_peak_caller.subgraphcollection import SubgraphCollection
from graph_peak_caller.areas import BinaryContinousAreas
from graph_peak_caller.peakscores import ScoredPeak
import offsetbasedgraph as obg
from graph_peak_caller.snarls import SnarlGraph, SnarlGraphBuilder
from graph_peak_caller.snarlmaps import LinearSnarlMap
import traceback
import warnings
import sys

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file, 'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = warn_with_traceback


def run_with_gam(gam_file_name, gam_control_file, vg_graph_file_name,
                 limit_to_chromosomes=False):
    retriever = SequenceRetriever.from_vg_graph("haplo1kg50-mhc.vg")
    logging.basicConfig(level=logging.INFO)
    logging.info("Running")

    # vg_graph = pyvg.Graph.create_from_file(vg_graph_file_name)
    # ob_graph = vg_graph.get_offset_based_graph()
    # ob_graph.to_file("obgraph")
    ob_graph = obg.GraphWithReversals.from_file("obgraph")
    #print(ob_graph.node_size(701))
    #return

    #builder = SnarlGraphBuilder.from_vg_snarls(ob_graph.copy(), "haplo1kg50-mhc.snarls")
    #snarlgraph = builder.build_snarl_graphs()
    # LinearSnarlMap(snarlgraph, ob_graph)
    #linear_map = LinearSnarlMap(snarlgraph, ob_graph)
    #linear_map.to_file("linear_map")
    linear_map = "linear_map"
    #snarlgraph._create_distance_dicts()

    #snarlgraph.to_file("haplo1kg50-mhc.snarlgraph")
    #snarlgraph = SnarlGraph.from_file("haplo1kg50-mhc.snarlgraph")

    #ob_graph = obg.GraphWithReversals.from_file("obgraph")
    graph_size = sum(block.length() for block in ob_graph.blocks.values())
    logging.info("Graph size: %d" % graph_size)
    logging.info("N nodes in graph: %d" % len(ob_graph.blocks))

    # print(ob_graph.blocks)
    reads_intervals = vg_gam_file_to_interval_collection(
         None, gam_file_name, ob_graph)

    control_intervals = vg_gam_file_to_interval_collection(
         None, gam_control_file, ob_graph)

    #experiment_info = callpeaks.ExperimentInfo(graph_size, 103, 50)
    experiment_info = callpeaks.ExperimentInfo(graph_size, 135, 36)
    caller = callpeaks.CallPeaks(
        ob_graph, reads_intervals, control_intervals,
        experiment_info=experiment_info,
        out_file_base_name="real_data_", has_control=True,
        linear_map=linear_map)
    caller.set_cutoff(0.10)
    caller.verbose = True
    caller.run()
    sequences = [retriever.get_interval_sequence(max_path)
                 for max_path in caller.max_paths]
    f = open("real_data_sequences", "w")
    i = 0
    for seq in sequences:
        f.write(">peak" + str(i) + "\n" + seq + "\n")
        i += 1


def run_from_max_paths_step(graph_file_name, pileup_file_name, read_length):
    graph = obg.Graph.from_file(graph_file_name)
    peaks = SparsePileup.from_bed_file(graph, pileup_file_name)
    peaks.fill_small_wholes(read_length)
    # final_track = peaks.remove_small_peaks(fragment_length)
    final_track = peaks
    peaks_as_subgraphs = final_track.to_subgraphs()
    p_values = SparsePileup.from_bed_file(graph, "real_data_q_values.bdg")
    # peaks_as_subgraphs = SubgraphCollection.from_file(graph, "real_data_peaks_as_subgraphs")
    binary_peaks = (BinaryContinousAreas.from_old_areas(peak) for peak in
                    peaks_as_subgraphs)
    scored_peaks = (ScoredPeak.from_peak_and_pileup(peak, p_values)
                    for peak in binary_peaks)
    max_paths = [scored_peak.get_max_path() for
                 scored_peak in scored_peaks]
    max_paths = [p for p in max_paths if p.length() > 136]
    IntervalCollection(max_paths).to_text_file(
                "real_data_max_paths")
    retriever = SequenceRetriever.from_vg_graph("haplo1kg50-mhc.vg")
    sequences = [retriever.get_interval_sequence(max_path)
                 for max_path in max_paths]
    f = open("real_data_sequences", "w")
    i = 0
    for seq in sequences:
        f.write(">peak" + str(i) + "\n" + seq + "\n")
        i += 1


def peak_sequences_to_fasta(vg_graph_file_name, peaks_file_name, fasta_file_name):

    sequence_retriever = SequenceRetriever.from_vg_graph(vg_graph_file_name)
    peaks = open(peaks_file_name)
    fasta_file = open(fasta_file_name, "w")
    for line in peaks.read():
        l = line.split()
        node = int(l[0])
        start = int(l[1])
        end = int(l[2])
        sequence = sequence_retriever.get_sequence_on_directed_node()




if __name__ == "__main__":
    dm_folder = "../graph_peak_caller/dm_test_data/"
    #run_from_max_paths_step("obgraph", "pre_postprocess.bed", 36)
    # run_from_max_paths_step("obgraph", "pre_postprocess.bed", 36)
    #cProfile.run('run_with_gam("ENCFF000WVQ_filtered.gam", "cactus-mhc.json")')
    #cProfile.run('run_with_gam("ENCFF001HNI_filtered_q60.gam", "ENCFF001HNS_filtered_q60.gam", "cactus-mhc.json")')
    #run_from_max_paths_step()
    #run_with_gam("ENCFF001HNI_filtered_q60.gam", "ENCFF001HNS_filtered_q60.gam", "cactus-mhc.json")
    #run_with_gam("ENCFF001HNI_filtered_q60.gam", "ENCFF001HNS_filtered_q60.gam", "haplo1kg50-mhc.json")
    run_with_gam("ENCFF001HNI_haplo1kg50-mhc_filtered_q30.gam", "ENCFF001HNS_haplo1kg50-mhc_filtered_q30.gam", "haplo1kg50-mhc.json")
