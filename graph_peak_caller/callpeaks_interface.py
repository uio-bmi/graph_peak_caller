import logging

import offsetbasedgraph as obg
from pyvg.conversion import vg_json_file_to_interval_collection
from pyvg.sequences import SequenceRetriever

from . import ExperimentInfo, Configuration, CallPeaks
from .multiplegraphscallpeaks import MultipleGraphsCallpeaks


def run_callpeaks(ob_graph,
                  sample_file_name, control_file_name,
                  vg_graph_file_name,
                  out_name="real_data_",
                  has_control=True,
                  limit_to_chromosomes=False,
                  fragment_length=135, read_length=36,
                  linear_map_file_name=False):

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
        skip_filter_duplicates=False, p_val_cutoff=0.05,
        graph_is_partially_ordered=True)

    caller = CallPeaks.run_from_intervals(
        ob_graph, sample_intervals, control_intervals,
        experiment_info=experiment_info,
        out_file_base_name=out_name, has_control=has_control,
        linear_map=linear_map_file_name,
        configuration=config
    )
    if vg_graph_file_name != "None":
        try:
            retriever = SequenceRetriever.from_vg_graph(vg_graph_file_name)
        except OSError:
            retriever = SequenceRetriever.from_vg_json_graph(vg_graph_file_name)

        caller.save_max_path_sequences_to_fasta_file("sequences.fasta", retriever)
    else:
        logging.info("Not saving max path sequences, since a vg graph/sequence retriever was not sent in")


def run_callpeaks_interface(args):
    logging.info("Read offset based graph")

    ob_graph = obg.GraphWithReversals.from_numpy_file(
        args.graph_file_name)

    run_callpeaks(
        ob_graph,
        args.sample_reads_file_name,
        args.control_reads_file_name,
        args.vg_graph_file_name,
        args.out_base_name,
        has_control=args.with_control == "True",
        fragment_length=int(args.fragment_length),
        read_length=int(args.read_length),
        linear_map_file_name=args.linear_map_base_name
    )


def run_callpeaks_whole_genome(args):
    logging.info("Running whole genome.")
    chromosomes = args.chromosomes.split(",")
    graph_file_names = [args.graphs_location + chrom for chrom in chromosomes]
    linear_map_file_names = [args.linear_maps_location + chrom
                             for chrom in chromosomes]
    vg_graphs = [args.vg_graphs_location + chrom + ".vg"
                 for chrom in chromosomes]
    sequence_retrievers = \
        (SequenceRetriever.from_vg_graph(fn) for fn in vg_graphs)

    if args.sample_reads_base_name.endswith(".intervalcollection"):
        sample_file_names = [
            args.sample_reads_base_name.replace("chrom", chrom)
            for chrom in chromosomes]
        control_file_names = [
            args.control_reads_base_name.replace("chrom", chrom)
            for chrom in chromosomes]
    else:
        sample_base_name = args.sample_reads_base_name.replace(".json", "_")
        control_base_name = args.control_reads_base_name.replace(".json", "_")
        sample_file_names = [sample_base_name + chrom + ".json"
                             for chrom in chromosomes]
        control_file_names = [control_base_name + chrom + ".json"
                              for chrom in chromosomes]

    fragment_length = int(args.fragment_length)
    genome_size = int(args.genome_size)
    min_background = int(args.unique_reads) * fragment_length / genome_size
    logging.info(
        "Computed min background signal to be %.3f using fragment length %f, "
        " %d unique reads, and genome size %d" % (int(min_background),
                                                  fragment_length,
                                                  int(args.unique_reads),
                                                  int(min_background)))

    caller = MultipleGraphsCallpeaks(
        chromosomes,
        graph_file_names,
        sample_file_names,
        control_file_names,
        linear_map_file_names,
        int(args.fragment_length),
        int(args.read_length),
        has_control=args.with_control == "True",
        sequence_retrievers=sequence_retrievers,
        out_base_name=args.out_base_name,
        stop_after_p_values=args.stop_after_p_values == "True",
        min_background=min_background
    )
    caller.run()


def run_callpeaks_whole_genome_from_p_values(args):
    logging.info("Running whole genome from p-values.")
    chromosome = args.chromosome
    chromosomes = [chromosome]
    graph_file_names = [args.graphs_location + chrom for chrom in chromosomes]
    vg_graphs = [args.graphs_location + chrom + ".vg" for chrom in chromosomes]
    sequence_retrievers = \
        (SequenceRetriever.from_vg_graph(fn) for fn in vg_graphs)

    caller = MultipleGraphsCallpeaks(
        chromosomes,
        graph_file_names,
        None,
        None,
        None,
        int(args.fragment_length),
        int(args.read_length),
        has_control=args.with_control == "True",
        sequence_retrievers=sequence_retrievers,
        out_base_name=args.out_base_name
    )
    caller.create_joined_q_value_mapping()
    caller.run_from_p_values(only_chromosome=chromosome)
