import pyvg.util
import pyvg.vg as vg
import offsetbasedgraph as obg
from pyvg.util import vg_gam_file_to_interval_collection, vg_json_file_to_interval_collection
from pyvg.sequences import SequenceRetriever
from graph_peak_caller.util import create_linear_map, create_ob_graph_from_vg
import logging
from graph_peak_caller.callpeaks import CallPeaks, ExperimentInfo, CallPeaksFromQvalues, Configuration
import argparse
import sys
from graph_peak_caller.sparsepileup import SparsePileup
import pickle
import subprocess
from pyvg.protoparser import json_file_to_obg_graph, json_file_to_obg_numpy_graph
from graph_peak_caller.peakcollection import Peak
import os
from collections import defaultdict
from memory_profiler import profile

logging.basicConfig(level=logging.INFO, format="%(asctime)s, %(levelname)s: %(message)s")


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
    #logging.info("N nodes in graph: %d" % len(ob_graph.blocks))

    experiment_info = ExperimentInfo(graph_size, fragment_length, read_length)
    config = Configuration(
        skip_read_validation=True, save_tmp_results_to_file=False,
        skip_filter_duplicates=True, p_val_cutoff=0.05,
        graph_is_partially_ordered=True)

    caller = CallPeaks.run_from_intervals(
        ob_graph, sample_intervals, control_intervals,
        experiment_info=experiment_info,
        out_file_base_name=out_name, has_control=has_control,
        linear_map=linear_map,
        configuration=config
    )
    retriever = SequenceRetriever.from_vg_graph(vg_graph_file_name)
    caller.save_max_path_sequences_to_fasta_file("sequences.fasta", retriever)


def run_with_gam(ob_graph,
                 gam_file_name, gam_control_file,
                 vg_graph_file_name,
                 out_name="real_data_",
                 has_control=True,
                 limit_to_chromosomes=False,
                 fragment_length=135, read_length=36,
                 linear_map_file_name = False):

    logging.info("Running from gam files")

    #ob_graph = obg.GraphWithReversals.from_file(ob_graph_file_name)
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

#@profile
def run_with_json(ob_graph,
                 json_file_name, json_control_file,
                 vg_graph_file_name,
                 out_name="real_data_",
                 has_control=True,
                 limit_to_chromosomes=False,
                 fragment_length=135, read_length=36,
                 linear_map_file_name = False):

    logging.info("Running from gam files")

    #ob_graph = obg.GraphWithReversals.from_file(ob_graph_file_name)
    reads_intervals = vg_json_file_to_interval_collection(
         None, json_file_name, ob_graph)

    control_intervals = vg_json_file_to_interval_collection(
         None, json_control_file, ob_graph)

    run_with_intervals(ob_graph, reads_intervals, control_intervals,
                       out_name=out_name, has_control=has_control,
                       vg_graph_file_name=vg_graph_file_name,
                       fragment_length=fragment_length,
                       read_length=read_length,
                       linear_map=linear_map_file_name)


def run_mhc_ctcf_example():
    create_ob_graph_from_vg("tests/mhc/graph.json", "tests/mhc/graph.obg")
    logging.info("Reading graph from file")
    ob_graph = obg.Graph.from_file("tests/mhc/graph.obg")
    create_linear_map(ob_graph, "tests/mhc/graph.snarls", "tests/mhc/linear_map.lm")

    run_with_gam(
        "tests/mhc/graph.obg",
        "tests/mhc/macs_remapped_mhc.gam",
        "tests/mhc/macs_remapped_mhc.gam",
        "tests/mhc/graph.vg",
        "tests/mhc/macs_reads_remapped_",
        has_control=False,
        fragment_length=135,
        read_length=36,
        linear_map_file_name="tests/mhc/linear_map.lm"
    )


def run_callpeaks_from_q_values(args):
    name = args.out_base_name
    logging.info("Reading obgraph from file")
    ob_graph = obg.GraphWithReversals.from_file(args.obg_file_name)
    logging.info("Number of nodes in graph: %d" %  len(ob_graph.blocks))

    logging.info("Creating q values pileup from file")
    q_values = SparsePileup.from_pickle(name + "q_values.pickle", ob_graph)
    logging.info("Getting touched nodes and exp info from file")
    with open(name + "touched_nodes.pickle", "rb") as f:
        touched_nodes = pickle.loads(f.read())

    experiment_info = ExperimentInfo.from_file(name + "experiment_info.pickle")
    logging.info("All prev data fetched")

    caller = CallPeaksFromQvalues(
        ob_graph,
        q_values,
        experiment_info,
        args.out_base_name,
        touched_nodes=touched_nodes,
        graph_is_partially_ordered=True
    )
    caller.callpeaks()


def intervals_to_fasta(args):
    logging.info("Getting sequence retriever")
    retriever = SequenceRetriever.from_vg_graph(args.vg_graph_file_name)
    logging.info("Getting intervals")
    intervals = obg.IntervalCollection.create_generator_from_file(args.intervals_file_name)
    logging.info("Writing to fasta")
    CallPeaksFromQvalues.intervals_to_fasta_file(intervals, args.out_file_name, retriever)


def run_callpeaks(args):
    logging.info("Creating offset based graph")
    out_name = args.out_base_name
    json_file_name = args.vg_json_graph_file_name
    obg_file_name = json_file_name.replace(".json", ".obg")

    n_nodes = 0  # args.n_nodes
    if not os.path.isfile(obg_file_name):
        ob_graph = json_file_to_obg_graph(json_file_name, n_nodes)
        logging.info("Writing ob graph to file")
        #ob_graph.to_numpy_files(obg_file_name)
        ob_graph.to_file(obg_file_name)
    else:
        logging.info("Reading graph from file (graph already existing on disk)")
        #ob_graph = obg.GraphWithReversals.from_numpy_files(obg_file_name)
        ob_graph = obg.GraphWithReversals.from_file(obg_file_name)


    if not os.path.isfile(args.linear_map_base_name + "_starts.pickle"):
        logging.info("Creating linear map")
        create_linear_map(ob_graph, args.vg_snarls_file_name, args.linear_map_base_name)
        logging.info("Linear map created")
    else:
        logging.info("Not creating linear map. Already existing")

    has_control = True
    if args.with_control == "False":
        has_control = False

    run_with_json(
        ob_graph,
        args.sample_reads_file_name,
        args.control_reads_file_name,
        args.vg_graph_file_name,
        args.out_base_name,
        has_control=has_control,
        fragment_length=int(args.fragment_length),
        read_length=int(args.read_length),
        linear_map_file_name=args.linear_map_base_name
    )


#@profile
def run_callpeaks_with_numpy_graph(args):
    logging.info("Read offset based graph")


    ob_graph = obg.GraphWithReversals.from_numpy_files(args.numpy_graph_base_name)

    if not os.path.isfile(args.linear_map_base_name + "_starts.pickle"):
        logging.info("Creating linear map")
        create_linear_map(ob_graph, args.vg_snarls_file_name, args.linear_map_base_name)
        logging.info("Linear map created")
    else:
        logging.info("Not creating linear map. Already existing")

    has_control = True
    if args.with_control == "False":
        has_control = False

    run_with_json(
        ob_graph,
        args.sample_reads_file_name,
        args.control_reads_file_name,
        args.vg_graph_file_name,
        args.out_base_name,
        has_control=has_control,
        fragment_length=int(args.fragment_length),
        read_length=int(args.read_length),
        linear_map_file_name=args.linear_map_base_name
    )


def linear_peaks_to_fasta(args):
    from benchmarking.nongraphpeaks import NonGraphPeakCollection
    collection = NonGraphPeakCollection.from_bed_file(args.linear_reads_file_name)
    #collection.filter_peaks_outside_region("chr6", 28510119, 33480577)
    collection.set_peak_sequences_using_fasta(fasta_file_location=args.fasta_file)
    collection.save_to_sorted_fasta(args.out_file_name)
    logging.info("Saved sequences to %s" % args.out_file_name)


def create_ob_graph(args):
    # Get node count
    #command = ["vg", "stats", "--node-count", args.vg_json_file_name]
    #result = subprocess.check_output(command, shell=True)
    #n_nodes = int(result)
    logging.info("Creating obgraph")
    ob_graph = json_file_to_obg_graph(args.vg_json_file_name, 0)
    logging.info("Writing ob graph to file")
    ob_graph.to_file(args.out_file_name)


def create_ob_numpy_graph(args):
    logging.info("Creating obgraph")
    ob_graph = json_file_to_obg_numpy_graph(args.vg_json_file_name, 0)

    logging.info("Writing ob graph to file")
    ob_graph.to_numpy_files(args.out_file_name)


def create_linear_map_interface(args):
    logging.info("Reading ob graph from file")
    ob_graph = obg.GraphWithReversals.from_file(args.obg_file_name)
    logging.info("Creating linear map")
    create_linear_map(ob_graph, args.vg_snarls_file_name, args.out_file_base_name)


def split_vg_json_reads_into_chromosomes(args):
    import json
    reads_base_name = args.vg_json_reads_file_name.split(".")[0]

    chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y"]
    chromosome_limits = {}
    logging.info("Found the following chromosome ranges:")
    for chrom in chromosomes:
        start_end = open(args.range_files_base_name + "node_range_" + chrom + ".txt")
        start_end= start_end.read().split(":")
        start = int(start_end[0])
        end = int(start_end[1])
        chromosome_limits[chrom] = (start, end)
        logging.info("   Chr%s: %d-%d" % (chrom, start, end))

    out_files = {chrom: open(reads_base_name + "_" + chrom + ".json", "w")
                 for chrom in chromosomes}

    reads_file = open(args.vg_json_reads_file_name)
    i = 0
    import re
    regex = re.compile(r"node_id\": ([0-9]+)")

    def get_mapped_chrom(node):
        for chrom in chromosomes:
            if node >= chromosome_limits[chrom][0] and node <= chromosome_limits[chrom][1]:
                mapped_chrom = chrom
                break
        assert mapped_chrom is not None, "Found no match for node id %d" % node
        return mapped_chrom

    for line in reads_file:
        if i % 100000 == 0:
            logging.info("Line #%d" % i)
        i += 1

        groups = regex.search(line).groups()
        if len(groups) > 0:
            node = int(groups[0])
            mapped_chrom = get_mapped_chrom(node)
            out_files[mapped_chrom].writelines([line])
        else:
            print("No groups fond")


    for file in out_files.values():
        file.close()

    logging.info("Done")


def concatenate_sequence_files(args):
    chromosomes = args.chromosomes.split(",")
    out_file_name = args.out_file_name

    all_fasta_entries = []
    for chromosome in chromosomes:
        logging.info("Processing chromosome %s" % chromosome)
        fasta_file = open("chr" + chromosome + "_sequences.fasta")
        for line in fasta_file:
            if line.startswith(">"):
                all_fasta_entries.append([line, None])
            else:
                # This is sequence, add to prev entry
                all_fasta_entries[-1][1] = line

    peaks = []
    for fasta_entry in all_fasta_entries:
        header = fasta_entry[0].rstrip()
        sequence = fasta_entry[1].rstrip()

        interval_json = "{%s" % header.split(" {")[1]
        interval = Peak.from_file_line(interval_json)
        interval.sequence = sequence
        peaks.append(interval)

    peaks = sorted(peaks, key=lambda s: -s.score)
    out_fasta = open(out_file_name, "w")
    i = 0
    for peak in peaks:
        out_fasta.writelines([">peak%d %s\n" % (i, peak.to_file_line())])
        out_fasta.writelines(["%s\n" % peak.sequence])
        i += 1
    out_fasta.close()

    logging.info("Wrote all peaks in sorted order to %s" % out_file_name)


def plot_motif_enrichment(args):
    from benchmarking.motifenrichment import MotifMatcher, plot_true_positives
    fasta1 = args.fasta1
    fasta2 = args.fasta2
    meme = args.meme_motif_file

    plot_true_positives(
        {
            "file1": fasta1,
            "file2": fasta2
        },
        meme,
        save_to_file=args.out_figure_file_name
    )


interface = \
{
    'callpeaks':
        {
            'help': 'Callpeaks',
            'arguments':
                [
                    ('vg_json_graph_file_name', "Json Graph file name (.json)"),
                    ('vg_graph_file_name', "Graph file name (.vg)"),
                    ('linear_map_base_name', "Set to desired base name. Will be used if exists, created if not."),
                    ('sample_reads_file_name', ' '),
                    ('control_reads_file_name', ' '),
                    ('with_control', 'True/False'),
                    ('out_base_name', 'eg experiment1_'),
                    ('fragment_length', ''),
                    ('read_length', '')
                ],
            'method': run_callpeaks
        },
    'callpeaks_with_numpy_graph':
        {
            'help': 'Callpeaks',
            'arguments':
                [
                    ('numpy_graph_base_name', ""),
                    ('vg_graph_file_name', "Graph file name (.vg)"),
                    ('linear_map_base_name', "Set to desired base name. Will be used if exists, created if not."),
                    ('sample_reads_file_name', ' '),
                    ('control_reads_file_name', ' '),
                    ('with_control', 'True/False'),
                    ('out_base_name', 'eg experiment1_'),
                    ('fragment_length', ''),
                    ('read_length', '')
                ],
            'method': run_callpeaks_with_numpy_graph
        },
    'callpeaks_from_qvalues':
        {
            'help': '...',
            'arguments':
                [
                    ('obg_file_name', 'Offsetbased graph file name'),
                    ('out_base_name', 'Out base name used on previous run (Used to fetch tmp files)')
                ],
            'method': run_callpeaks_from_q_values
        },
    'intervals_to_fasta':
        {
            'help': 'Get sequences for intervals in interval file. Write to fasta',
            'arguments':
                [
                    ('vg_graph_file_name', ''),
                    ('intervals_file_name', ''),
                    ('out_file_name', '')
                ],
            'method': intervals_to_fasta
        },
    'linear_peaks_to_fasta':
        {
            'help': 'Converts a linear peaks file (eg. from macs) to fasta',
            'arguments':
                [
                    ('linear_reads_file_name', ''),
                    ('fasta_file', ''),
                    ('out_file_name', '')
                ],
            'method': linear_peaks_to_fasta
        },
    'create_ob_graph':
        {
            'help': 'Creates and stores an obgraph from a vg json graph',
            'arguments':
                [
                    ('vg_json_file_name', 'Vg json file name (created by running vg view -Vj graph.vg > graph.json'),
                    ('out_file_name', 'E.g. graph.obg')
                ],
            'method': create_ob_graph
        },
    'create_ob_numpy_graph':
        {
            'help': 'Creates and stores an obgraph from a vg json graph (using numpy data structures)',
            'arguments':
                [
                    ('vg_json_file_name', 'Vg json file name (created by running vg view -Vj graph.vg > graph.json'),
                    ('out_file_name', 'E.g. graph.obg')
                ],
            'method': create_ob_numpy_graph
        },
    'create_linear_map':
        {
            'help': 'Creates a linear map using a vg snarls file and an ob graph',
            'arguments':
                [
                    ('obg_file_name', ''),
                    ('vg_snarls_file_name', ''),
                    ('out_file_base_name', ''),
                ],
            'method': create_linear_map_interface
        },
    'split_vg_json_reads_into_chromosomes':
        {
            'help': "Splits intervals from interval collection into one file for each chromsome."
                    "Requires node_range_[chrom ID].txt to exist for each chromosome.",
            'arguments':
                [
                    ('vg_json_reads_file_name', ''),
                    ('range_files_base_name', 'Base name, e.g. dir, to range files')
                ],
            'method': split_vg_json_reads_into_chromosomes
        },
    'concatenate_sequence_files':
        {
            'help': "Merge multiple *_sequence.fasta files from the peak caller into one single sorted file.",
            'arguments':
                [
                    ('chromosomes', 'comma delimted, e.g 1,2,3, used to fetch files of type chr1_sequences.fasta, ...'),
                    ('out_file_name', '')
                ],
            'method': concatenate_sequence_files
        },
    'plot_motif_enrichment':
        {
            'help': "Plots motif enrichments using fimo. Requires fimo to be installed and in path.",
            'arguments':
                [
                    ('fasta1', ''),
                    ('fasta2', ''),
                    ('meme_motif_file', 'something.meme'),
                    ('out_figure_file_name', '')
                ],
            'method': plot_motif_enrichment
        }
}

# Create parser
parser = argparse.ArgumentParser(
    description='Graph peak caller')
subparsers = parser.add_subparsers(help='Subcommands')

for command in interface:
    example = ""
    if "example_run" in interface[command]:
        example = "\nExample: " + interface[command]["example_run"]

    subparser = subparsers.add_parser(command,
                            help=interface[command]["help"] + example)
    for argument, help in interface[command]["arguments"]:
        subparser.add_argument(argument, help=help)
    subparser.set_defaults(func=interface[command]["method"])

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()
if hasattr(args, 'func'):
    args.func(args)
else:
    parser.help()


"""
Run on server from start:
python3 ../../dev/graph_peak_caller/graph_peak_caller.py callpeaks graph.json graph.vg graph.snarls filtered_short.gam filtered_short.gam False run1/ 135 36 23739138

Using chr1, run4 and full dataset, without control:
cd data/chr1
python3 ../../dev/graph_peak_caller/graph_peak_caller.py callpeaks graph.json graph.vg graph.snarls ctcf_filtered_r0.97_e0.3-2.gam ctcf_filtered_r0.97_e0.3-2.gam False run4/ 135 36 23739138


Chr 19:
cd data/chr19
python3 ../../dev/graph_peak_caller/graph_peak_caller.py callpeaks graph.json graph.vg graph.snarls filtered_short.gam filtered_short.gam False run1/ 135 36 6373453

"""
"""
Lrc_kir local:
python3 ../../graph_peak_caller.py callpeaks graph.json graph.vg linear_map reads.json reads.json False test_ 136 35

LRC med numpy graph:
python3 ../../graph_peak_caller.py callpeaks_with_numpy_graph graph graph.vg linear_map reads.json reads.json False test_ 136 35


python3 ../../dev/graph_peak_caller/graph_peak_caller.py callpeaks graph.json graph.vg graph.snarls filtered.gam filtered.gam False run1/ 135 36 23739138



python3 ../../graph_peak_caller.py callpeaks_from_qvalues graph.obg run1/q_values.bdg 136 35 1098808 run1/

Run on server from q values:
python3 ../../dev/graph_peak_caller/graph_peak_caller.py callpeaks_from_qvalues graph.obg run2/

Run lrc_kir from qvalues:
/../graph_peak_caller.py callpeaks_from_qvalues graph.obg test2_


python3 ../graph_peak_caller.py linear_peaks_to_fasta macs_with_control_peaks_chr1.narrowPeak macs_with_control_sequences_chr1.fasta



Max TF chr1 på server
cd ~/data/chr1
python3 ../../dev/graph_peak_caller/graph_peak_caller.py callpeaks graph.json graph.vg graph.snarls ~/data/tfs/max/filtered_r1.0.gam ~/data/tfs/max/filtered_r1.0.gam False max_without_control/ 183 50 23739138

"""
