import logging

import offsetbasedgraph as obg
from pyvg.conversion import vg_json_file_to_interval_collection
from pyvg.sequences import SequenceRetriever

from . import ExperimentInfo, Configuration, CallPeaks
from .multiplegraphscallpeaks import MultipleGraphsCallpeaks
import os
from graph_peak_caller.util import create_linear_map

def run_callpeaks(ob_graph,
                  sample_file_name, control_file_name,
                  sequence_graph,
                  out_name="real_data_",
                  has_control=True,
                  limit_to_chromosomes=False,
                  fragment_length=135, read_length=36,
                  linear_map_file_name=False,
                  qval_threshold=0.05):

    logging.info("Running from gam files")

    # Detect file format
    if isinstance(sample_file_name, obg.IntervalCollection):
        sample_intervals = sample_file_name
        control_intervals = control_file_name
    elif sample_file_name.endswith(".json"):
        sample_intervals = vg_json_file_to_interval_collection(
            sample_file_name, ob_graph)
        control_intervals = vg_json_file_to_interval_collection(
            control_file_name, ob_graph)
    else:
        try:
            sample_intervals = obg.IntervalCollection.from_file(
                sample_file_name, graph=ob_graph)
            control_intervals = obg.IntervalCollection.from_file(
                control_file_name, graph=ob_graph)
        except OSError:
            sample_intervals = obg.IntervalCollection.from_file(
                sample_file_name, graph=ob_graph, text_file=True)
            control_intervals = obg.IntervalCollection.from_file(
                control_file_name, graph=ob_graph, text_file=True)

    graph_size = ob_graph.number_of_basepairs()
    logging.info("Number of base pairs in graph: %d" % graph_size)

    experiment_info = ExperimentInfo(graph_size, fragment_length, read_length)
    config = Configuration(
        skip_read_validation=True, save_tmp_results_to_file=True,
        skip_filter_duplicates=False, p_val_cutoff=qval_threshold,
        graph_is_partially_ordered=True)

    caller = CallPeaks.run_from_intervals(
        ob_graph, sample_intervals, control_intervals,
        experiment_info=experiment_info,
        out_file_base_name=out_name, has_control=has_control,
        linear_map=linear_map_file_name,
        configuration=config
    )
    if sequence_graph != None:
        caller.save_max_path_sequences_to_fasta_file("sequences.fasta", sequence_graph)
    else:
        logging.info("Not saving max path sequences, since a sequence graph was not found.")


def find_or_create_linear_map(graph, linear_map_name):

    if os.path.isfile(linear_map_name):
        logging.info("Found linear map %s. "
                     "Will be used in peak calling." % linear_map_name)
    else:
        logging.info("No linear map provided, and none found. Will create now.")
        create_linear_map(graph, linear_map_name)
        logging.info("Done creating linear map")


def run_callpeaks_interface(args):
    logging.info("Read offset based graph")

    ob_graph = args.graph
    control = args.sample
    has_control = False
    if args.control is not None:
        control = args.control
        has_control = True

    out_name = args.out_name if args.out_name is not None else ""
    linear_map_name = args.linear_map
    if linear_map_name is None:
        linear_map_name = args.graph_file_name.split(".")[0] + "_linear_map.npz"
        find_or_create_linear_map(ob_graph, linear_map_name)

    qval = 0.05 if args.q_value_threshold is None else float(args.q_value_thresholds)
    run_callpeaks(
        ob_graph,
        args.sample,
        control,
        args.sequence_graph,
        out_name,
        has_control=has_control,
        fragment_length=int(args.fragment_length),
        read_length=int(args.read_length),
        linear_map_file_name=linear_map_name,
        qval_threshold=qval
    )


def run_callpeaks_whole_genome(args):
    logging.info("Running whole genome.")
    chromosomes = args.chromosomes.split(",")
    graph_file_names = [args.data_dir + "/" + chrom for chrom in chromosomes]

    linear_map_file_names = []
    for i, chrom in enumerate(chromosomes):
        linear_map_name = args.data_dir + "/" + chrom + "_linear_map.npz"
        if not os.path.isfile(linear_map_name):
            logging.warning("Did not find linear map for "
                            "chromosome %s. Will create." % chrom)
            create_linear_map(obg.Graph.from_numpy_file(graph_file_names[i]), linear_map_name)
        else:
            logging.info("Found linear map %s that will be used." % linear_map_name)
        linear_map_file_names.append(linear_map_name)

    sequence_retrievers = \
        (obg.SequenceGraph.from_file(args.data_dir + "/" + chrom + ".nobg.sequences")
         for chrom in chromosomes)

    if args.sample.endswith(".intervalcollection"):
        sample_file_names = [args.sample.replace("chrom", chrom) for chrom in chromosomes]

    else:
        sample_base_name = args.sample.replace(".json", "_")
        sample_file_names = [sample_base_name + chrom + ".json" for chrom in chromosomes]

    if args.control is not None:
        if args.control.endswith(".intervalcollection"):
            control_file_names = [
                args.control_reads_base_name.replace("chrom", chrom)
                for chrom in chromosomes]
        else:
            control_base_name = args.control.replace(".json", "_")
            control_file_names = [control_base_name + chrom + ".json"
                                  for chrom in chromosomes]
    else:
        control_file_names = sample_file_names.copy()

    fragment_length = int(args.fragment_length)
    genome_size = int(args.genome_size)
    min_background = int(args.unique_reads) * fragment_length / genome_size
    logging.info(
        "Computed min background signal to be %.3f using fragment length %f, "
        " %d unique reads, and genome size %d" % (min_background,
                                                  fragment_length,
                                                  int(args.unique_reads),
                                                  int(genome_size)))

    out_name = args.out_name if args.out_name is not None else ""

    caller = MultipleGraphsCallpeaks(
        chromosomes,
        graph_file_names,
        sample_file_names,
        control_file_names,
        linear_map_file_names,
        int(args.fragment_length),
        int(args.read_length),
        has_control=args.control is not None,
        sequence_retrievers=sequence_retrievers,
        out_base_name=out_name,
        stop_after_p_values=args.stop_after_p_values == "True",
        min_background=min_background
    )
    caller.run()


def run_callpeaks_whole_genome_from_p_values(args):
    logging.info("Running whole genome from p-values.")
    chromosome = args.chromosome
    chromosomes = [chromosome]
    graph_file_names = [args.data_dir + chrom for chrom in chromosomes]
    sequence_retrievers = \
        (obg.SequenceGraph.from_file(args.data_dir + "/" + chrom + ".nobg.sequences") for chrom in chromosomes)
    out_name = args.out_name if args.out_name is not None else ""

    caller = MultipleGraphsCallpeaks(
        chromosomes,
        graph_file_names,
        None,
        None,
        None,
        int(args.fragment_length),
        int(args.read_length),
        has_control=None,
        sequence_retrievers=sequence_retrievers,
        out_base_name=out_name
    )
    caller.create_joined_q_value_mapping()
    caller.run_from_p_values(only_chromosome=chromosome)
