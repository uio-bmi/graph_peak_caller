import pyvg.util
import pyvg.vg as vg
import offsetbasedgraph as obg
from pyvg.util import vg_gam_file_to_interval_collection
from pyvg.sequences import SequenceRetriever
from graph_peak_caller.util import create_linear_map, create_ob_graph_from_vg
import logging
from graph_peak_caller.callpeaks import CallPeaks, ExperimentInfo
import argparse
import sys

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
    logging.info("N nodes in graph: %d" % len(ob_graph.blocks))

    experiment_info = ExperimentInfo(graph_size, fragment_length, read_length)
    caller = CallPeaks(
        ob_graph, sample_intervals, control_intervals,
        experiment_info=experiment_info,
        out_file_base_name=out_name, has_control=has_control,
        linear_map=linear_map)
    caller.set_cutoff(0.05)
    caller.verbose = True
    caller.run()
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


def run_callpeaks(args):
    logging.info("Creating offset based graph")
    from pyvg.protoparser import json_file_to_obg_graph
    import os
    out_name = args.out_base_name
    json_file_name = args.vg_json_graph_file_name
    obg_file_name = json_file_name.replace(".json", ".obg")

    if not os.path.isfile(obg_file_name + ".npy"):
        ob_graph = json_file_to_obg_graph(json_file_name, int(args.n_nodes))
        logging.info("Writing ob graph to file")
        ob_graph.to_numpy_files(obg_file_name)
    else:
        logging.info("Reading graph from file (graph already existing on disk")
        ob_graph = obg.GraphWithReversals.from_numpy_files(obg_file_name)

    if not os.path.isfile(out_name + "linear_map" + "starts.pickle"):
        logging.info("Creating linear map")
        create_linear_map(ob_graph, args.vg_snarls_file_name, out_name + "linear_map")
        logging.info("Linear map created")
    else:
        logging.info("Not creating linear map. Already existing")

    has_control = True
    if args.with_control == "False":
        has_control = False

    run_with_gam(
        ob_graph,
        args.sample_reads_file_name,
        args.control_reads_file_name,
        args.vg_graph_file_name,
        args.out_base_name,
        has_control=has_control,
        fragment_length=int(args.fragment_length),
        read_length=int(args.read_length),
        linear_map_file_name=out_name + "linear_map"
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
                    ('vg_snarls_file_name', "Snarls file name"),
                    ('sample_reads_file_name', ' '),
                    ('control_reads_file_name', ' '),
                    ('with_control', 'True/False'),
                    ('out_base_name', 'eg experiment1_'),
                    ('fragment_length', ''),
                    ('read_length', ''),
                    ('n_nodes', 'Number of nodes in grap')
                ],
            'method': run_callpeaks
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
python3 ../../dev/graph_peak_caller/graph_peak_caller.py callpeaks graph.json graph.vg graph.snarls filtered_short.gam filtered_short.gam False run1/ 135 36 23739138

Chr 19:
python3 ../../dev/graph_peak_caller/graph_peak_caller.py callpeaks graph.json graph.vg graph.snarls filtered_short.gam filtered_short.gam False run1/ 135 36 6373453

"""
"""
python3 ../../graph_peak_caller.py callpeaks graph.json graph.vg graph.snarls ctcf_filtered_r1.0.2.gam ctcf_filtered_r1.0.2.gam False test_ 136 35 112342

python3 ../../dev/graph_peak_caller/graph_peak_caller.py callpeaks graph.json graph.vg graph.snarls filtered.gam filtered.gam False run1/ 135 36 23739138


"""
