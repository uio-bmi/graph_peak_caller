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
from graph_peak_caller.peakcollection import PeakCollection
import traceback
import warnings
import sys


def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file, 'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = warn_with_traceback


def create_linear_map(ob_graph):
    builder = SnarlGraphBuilder.from_vg_snarls(ob_graph.copy(), "haplo1kg50-mhc.snarls")
    snarlgraph = builder.build_snarl_graphs()
    LinearSnarlMap(snarlgraph, ob_graph)
    linear_map = LinearSnarlMap(snarlgraph, ob_graph)
    linear_map.to_file("linear_map")

    # vg_graph = pyvg.Graph.create_from_file(vg_graph_file_name)
    # ob_graph = vg_graph.get_offset_based_graph()
    # ob_graph.to_file("obgraph")
    ob_graph = obg.GraphWithReversals.from_file("graph.obg")
    #print(ob_graph.node_size(701))
    #return

    builder = SnarlGraphBuilder.from_vg_snarls(
        ob_graph.copy(),
        "haplo1kg50-mhc.snarls")
    snarlgraph = builder.build_snarl_graphs()
    # LinearSnarlMap(snarlgraph, ob_graph)
    linear_map = LinearSnarlMap(snarlgraph, ob_graph)
    linear_map.to_file("haplo1kg50-mhc.lm")
    linear_map = "haplo1kg50-mhc.lm"
    #snarlgraph._create_distance_dicts()


def run_with_intervals(sample_intervals, control_intervals, out_name, has_control=True):
    logging.info("Running from intervals")
    retriever = SequenceRetriever.from_vg_graph("haplo1kg50-mhc.vg")
    ob_graph = obg.GraphWithReversals.from_file("graph.obg")
    graph_size = sum(block.length() for block in ob_graph.blocks.values())
    logging.info("Graph size: %d" % graph_size)
    logging.info("N nodes in graph: %d" % len(ob_graph.blocks))

    linear_map = "haplo1kg50-mhc.lm"
    experiment_info = callpeaks.ExperimentInfo(graph_size, 135, 36)
    caller = callpeaks.CallPeaksWRawReads(
        ob_graph, sample_intervals, control_intervals,
        experiment_info=experiment_info,
        out_file_base_name=out_name, has_control=has_control,
        linear_map=linear_map)
    caller.set_cutoff(0.05)
    caller.verbose = True
    caller.run()
    retriever = SequenceRetriever.from_vg_graph("haplo1kg50-mhc.vg")
    caller.save_max_path_sequences_to_fasta_file("sequences.fasta", retriever)


def run_with_gam(gam_file_name, gam_control_file, vg_graph_file_name,
                 out_name="real_data_",
                 has_control=True,
                 limit_to_chromosomes=False):
    logging.basicConfig(level=logging.INFO)
    logging.info("Running from gam files")

    ob_graph = obg.GraphWithReversals.from_file("graph.obg")
    # print(ob_graph.blocks)
    reads_intervals = vg_gam_file_to_interval_collection(
         None, gam_file_name, ob_graph)

    control_intervals = vg_gam_file_to_interval_collection(
         None, gam_control_file, ob_graph)

    run_with_intervals(reads_intervals, control_intervals,
                       out_name=out_name, has_control=has_control)


def run_from_max_paths_step(graph_file_name, pileup_file_name, raw_pileup_file_name):
    ob_graph = obg.GraphWithReversals.from_file("obgraph")
    graph_size = sum(block.length() for block in ob_graph.blocks.values())
    experiment_info = callpeaks.ExperimentInfo(graph_size, 135, 36)
    q_values = SparsePileup.from_bed_graph(ob_graph, pileup_file_name)
    raw_pileup = SparsePileup.from_bed_graph(ob_graph, raw_pileup_file_name)
    fromqvalues = callpeaks.CallPeaksFromQvalues(
        ob_graph, q_values, experiment_info, "laststep", raw_pileup=raw_pileup)

    fromqvalues.callpeaks()
    retriever = SequenceRetriever.from_vg_graph("haplo1kg50-mhc.vg")
    caller.save_max_path_sequences_to_fasta_file("sequences.fasta", retriever)
    # fragment_length = 135
    # graph = obg.Graph.from_file(graph_file_name)
    # peaks = SparsePileup.from_bed_file(graph, pileup_file_name)
    # peaks.fill_small_wholes(read_length)
    # final_track = peaks.remove_small_peaks(fragment_length)
    # peaks_as_subgraphs = final_track.to_subgraphs()
    # peaks_as_subgraphs.to_file(
    #     "last_step_" + "peaks_as_subgraphs")
    # 
    # p_values = SparsePileup.from_bed_file(graph, "real_data_q_values.bdg")
    # binary_peaks = (BinaryContinousAreas.from_old_areas(peak) for peak in
    #                 peaks_as_subgraphs)
    # scored_peaks = (ScoredPeak.from_peak_and_pileup(peak, p_values)
    #                 for peak in binary_peaks)
    # max_paths = [scored_peak.get_max_path() for
    #              scored_peak in scored_peaks]
    # max_paths.sort(key=lambda p: p.score, reverse=True)
    # PeakCollection(max_paths).to_file(
    #     "last_step_" + "max_paths", text_file=True)
    # 
    # # IntervalCollection(max_paths).to_text_file(
    # #             "last_step_max_paths")


def get_sequences(path_file):
    max_paths = PeakCollection.from_file(path_file, True)
    retriever = SequenceRetriever.from_vg_graph("haplo1kg50-mhc.vg")
    sequences = [retriever.get_interval_sequence(max_path)
                 for max_path in max_paths]
    f = open("tmp_sequences", "w")
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


def create_ob_graph_from_vg(vg_json_graph_file_name, ob_graph_file_name="graph.obg"):
    vg_graph = pyvg.Graph.create_from_file(vg_json_graph_file_name)
    ob_graph = vg_graph.get_offset_based_graph()
    ob_graph.to_file(ob_graph_file_name)
    logging.info("Wrote obgraph to %s" % ob_graph_file_name)


def run_ctcf_example():
     run_with_gam("ENCFF001HNI_haplo1kg50-mhc_filtered_q50.gam",
                 "ENCFF001HNI_haplo1kg50-mhc_filtered_q50.gam",
                 "haplo1kg50-mhc.json",
                 out_name="hack_",
                  has_control=False)


def run_ctcf_example_w_control():
     run_with_gam("ENCFF001HNI_haplo1kg50-mhc_filtered_q50.gam",
                  "ENCFF001HNS_haplo1kg50-mhc_filtered_q50.gam",
                  "haplo1kg50-mhc.json",
                  out_name="ctcf_q50_with_control_",
                  has_control=True)

if __name__ == "__main__":
    # get_sequences("laststepmax_paths.intervalcollection")
    run_ctcf_example_w_control()
    exit()
    # run_from_max_paths_step("obgraph", "ctcf_q50_with_control_q_values.bdg",
    # "ctcf_q50_with_control_raw_track.bdg")
    exit()

    dm_folder = "../graph_peak_caller/dm_test_data/"
    # ob_graph = obg.GraphWithReversals.from_file("obgraph")
    # create_linear_map(ob_graph)
    # create_ob_graph_from_vg("haplo1kg50-mhc.json")
    ob_graph = obg.GraphWithReversals.from_file("graph.obg")
    #create_linear_map(ob_graph)
    
    #cProfile.run('run_with_gam("ENCFF000WVQ_filtered.gam", "cactus-mhc.json")')
    #cProfile.run('run_with_gam("ENCFF001HNI_filtered_q60.gam", "ENCFF001HNS_filtered_q60.gam", "cactus-mhc.json")')
    #run_from_max_paths_step()
    #run_with_gam("ENCFF001HNI_filtered_q60.gam", "ENCFF001HNS_filtered_q60.gam", "cactus-mhc.json")
    #run_with_gam("ENCFF001HNI_filtered_q60.gam", "ENCFF001HNS_filtered_q60.gam", "haplo1kg50-mhc.json")

    #run_with_intervals(
    #    sample_intervals=IntervalCollection.from_file("sample_linear_reads.intervals", graph=ob_graph),
    #    control_intervals=IntervalCollection.from_file("control_linear_reads.intervals", graph=ob_graph)
    #)

    #run_with_gam("ENCFF001HNI_haplo1kg50-mhc_filtered_q30.gam", "ENCFF001HNS_haplo1kg50-mhc_filtered_q30.gam", "haplo1kg50-mhc.json")
    # run_with_gam("ctcf_mhc.gam", "ctcf_control_mhc.gam", "haplo1kg50-mhc.json")

