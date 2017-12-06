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
from graph_peak_caller.util import create_linear_map, create_ob_graph_from_vg
import traceback
import warnings
import sys


def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file, 'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = warn_with_traceback


def run_with_intervals(ob_graph,
                       sample_intervals,
                       control_intervals,
                       out_name,
                       has_control=True,
                       vg_graph_file_name="haplo1kg50-mhc.vg",
                       fragment_length=135,
                       read_length=36,
                       linear_map="haplo1kg50-mhc.lm"):
    logging.info("Running from intervals")
    graph_size = sum(block.length() for block in ob_graph.blocks.values())
    logging.info("Graph size: %d" % graph_size)
    logging.info("N nodes in graph: %d" % len(ob_graph.blocks))

    experiment_info = callpeaks.ExperimentInfo(graph_size, fragment_length, read_length)
    caller = callpeaks.CallPeaks(
        ob_graph, sample_intervals, control_intervals,
        experiment_info=experiment_info,
        out_file_base_name=out_name, has_control=has_control,
        linear_map=linear_map)
    caller.set_cutoff(0.05)
    caller.verbose = True
    caller.run()
    retriever = SequenceRetriever.from_vg_graph(vg_graph_file_name)
    caller.save_max_path_sequences_to_fasta_file("sequences.fasta", retriever)


def run_with_gam(ob_graph_file_name,
                 gam_file_name, gam_control_file,
                 vg_graph_file_name,
                 out_name="real_data_",
                 has_control=True,
                 limit_to_chromosomes=False,
                 fragment_length=135, read_length=36,
                 linear_map_file_name = False):

    logging.basicConfig(level=logging.INFO)
    logging.info("Running from gam files")

    ob_graph = obg.GraphWithReversals.from_file(ob_graph_file_name)
    reads_intervals = vg_gam_file_to_interval_collection(
         None, gam_file_name, ob_graph)

    control_intervals = vg_gam_file_to_interval_collection(
         None, gam_control_file, ob_graph)

    run_with_intervals(ob_graph, reads_intervals, control_intervals,
                       out_name=out_name, has_control=has_control,
                       vg_graph_file_name=vg_graph_file_name,
                       fragment_length=fragment_length,
                       read_length=read_length,
                       linear_map=linear_map_file_name)


def run_from_max_paths_step(graph_file_name, pileup_file_name, raw_pileup_file_name=None):
    raise Exception("Deprecated. Needs to be updated")
    ob_graph = obg.GraphWithReversals.from_file("obgraph")
    graph_size = sum(block.length() for block in ob_graph.blocks.values())
    experiment_info = callpeaks.ExperimentInfo(graph_size, 135, 36)
    q_values = SparsePileup.from_bed_graph(ob_graph, pileup_file_name)
    if raw_pileup_file_name is not None:
        raw_pileup = SparsePileup.from_bed_graph(ob_graph, raw_pileup_file_name)
    else:
        raw_pileup = None
    fromqvalues = callpeaks.CallPeaksFromQvalues(
        ob_graph, q_values, experiment_info, "laststep_", raw_pileup=raw_pileup, cutoff=0.1)

    fromqvalues.callpeaks()
    retriever = SequenceRetriever.from_vg_graph("haplo1kg50-mhc.vg")
    fromqvalues.save_max_path_sequences_to_fasta_file("sequences.fasta",
                                                      retriever)



def run_with_macs_filtered_reads():
    run_with_gam("ctcf_mhc.gam",
                 "ctcf_mhc.gam",
                 "haplo1kg50-mhc.json",
                 out_name="ctcf_macs_reads_without_control_",
                  has_control=False)


def run_with_reads_filtered_outside():
    run_with_gam("ctcf_mhc_data/ctcf_without_outside2_filtered_q30.gam",
                 "ctcf_mhc_data/ctcf_without_outside2_filtered_q30.gam",
                 "haplo1kg50-mhc.json",
                 out_name="ctcf_filtered_outside_",
                  has_control=False)


def run_with_linear_reads_moved_to_graph_without_control():
    ob_graph = obg.GraphWithReversals.from_file("graph.obg")
    run_with_intervals(
        sample_intervals=IntervalCollection.from_file("sample_linear_reads.intervals", graph=ob_graph),
        control_intervals=IntervalCollection.from_file("sample_linear_reads.intervals", graph=ob_graph),
        has_control=False,
        out_name="linear_reads_moved_to_graph_"
    )


def run_macs_reads_remapped_without_control():
    run_with_gam("ctcf_mhc.gam",
                 "ctcf_mhc.gam",
                 "haplo1kg50-mhc.json",
                 out_name="ctcf_macs_reads_remapped_without_control_",
                  has_control=False)


def run_with_macs_filtered_reads_w_control():
    run_with_gam("ctcf_mhc.gam",
                 "ctcf_control_mhc.gam",
                 "haplo1kg50-mhc.json",
                 out_name="ctcf_macs_reads_with_control_",
                  has_control=False)


def run_ctcf_example():
    file_name = "ENCFF001HNI_haplo1kg50-mhc_filtered_q50.gam"
    # file_name = "ENCFF001HNI_filtered_q60_r099.gam",
    # file_name = "ENCFF001HNI_filtered_q60_r30.gam"
    run_with_gam(file_name,
                 file_name,
                 "haplo1kg50-mhc.json",
                 out_name="ctcf_q60_subsampled_",
                 has_control=False)


def run_ctcf_example_w_control():
    run_with_gam("ENCFF001HNI_haplo1kg50-mhc_filtered_q50.gam",
                 "ENCFF001HNS_haplo1kg50-mhc_filtered_q50.gam",
                 "haplo1kg50-mhc.json",
                 out_name="ctcf_q50_with_control_2_",
                 has_control=True)

def run_srf_example():
    run_with_gam("vgdata/srf_filtered_r1.0_2.gam",
                 "vgdata/srf_filtered_r1.0_2.gam",
                 "haplo1kg50-mhc.json",
                 out_name="srf_",
                 has_control=False,
                 fragment_length=161,
                 read_length=50)

def run_lrc_kir_ctcf_example():
    #create_ob_graph_from_vg("lrc_kir/graph.json", "lrc_kir/graph.obg")
    #ob_graph = obg.Graph.from_file("lrc_kir/graph.obg")
    #create_linear_map(ob_graph, "lrc_kir/graph.snarls", "lrc_kir/lrc_kir.lm")

    run_with_gam(
        "lrc_kir/graph.obg",
        "lrc_kir/macs_reads_remapped.gam",
        "lrc_kir/macs_reads_remapped.gam",
        "lrc_kir/graph.vg",
        "lrc_kir/macs_reads_remapped_",
        has_control=False,
        fragment_length=135,
        read_length=36,
        linear_map_file_name="lrc_kir/lrc_kir.lm"
    )

def run_mhc_ctcf_example():
    create_ob_graph_from_vg("mhc/graph.json", "mhc/graph.obg")
    ob_graph = obg.Graph.from_file("mhc/graph.obg")
    create_linear_map(ob_graph, "mhc/graph.snarls", "mhc/linear_map.lm")

    run_with_gam(
        "mhc/graph.obg",
        "mhc/macs_remapped_mhc.gam",
        "mhc/macs_remapped_mhc.gam",
        "mhc/graph.vg",
        "mhc/macs_reads_remapped_",
        has_control=False,
        fragment_length=135,
        read_length=36,
        linear_map_file_name="mhc/linear_map.lm"
    )

def run_lrc_kir_using_macs_reads():
    ob_graph = obg.Graph.from_file("lrc_kir/graph.obg")
    reads1 = IntervalCollection.from_file("lrc_kir/macs_reads_on_graph.intervals", text_file=True, graph=ob_graph)
    reads2 = IntervalCollection.from_file("lrc_kir/macs_reads_on_graph.intervals", text_file=True, graph=ob_graph)

    run_with_intervals(ob_graph,
                       reads1, reads2,
                       "lrc_kir_using_macs_reads_",
                       has_control=False,
                       vg_graph_file_name="lrc_kir/graph.vg",
                       fragment_length=135,
                       read_length=36,
                       linear_map="lrc_kir/lrc_kir.lm"
    )

def run_mhc_using_macs_reads():
    ob_graph = obg.Graph.from_file("mhc/graph.obg")
    reads1 = IntervalCollection.from_file("mhc/macs_reads_on_graph.intervals", text_file=True, graph=ob_graph)
    reads2 = IntervalCollection.from_file("mhc/macs_reads_on_graph.intervals", text_file=True, graph=ob_graph)

    run_with_intervals(ob_graph,
                       reads1, reads2,
                       "mhc_using_macs_reads_",
                       has_control=False,
                       vg_graph_file_name="mhc/graph.vg",
                       fragment_length=135,
                       read_length=36,
                       linear_map="mhc/linear_map.lm"
    )

def run_chr1_example():
    create_ob_graph_from_vg("vgdata/chr1/graph.json", "vgdata/chr1/graph.obg")
    ob_graph = obg.Graph.from_file("vgdata/chr1/graph.obg")
    create_linear_map(ob_graph, "vgdata/chr1/graph.snarls", "vgdata/chr1/linear_map.lm")


def run_from_q_values(out_name):
    pileup_name = out_name + "q_values.bdg"
    run_from_max_paths_step("graph.obg", pileup_name)

if __name__ == "__main__":
    run_chr1_example()
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

