import logging
import os

import offsetbasedgraph as obg
from pyvg.conversion import vg_json_file_to_interval_collection

from . import Configuration, CallPeaks
from .multiplegraphscallpeaks import MultipleGraphsCallpeaks
from .util import create_linear_map
from .peakfasta import PeakFasta
from .reporter import Reporter
from .intervals import UniqueIntervals
import sys


def get_confiugration(args):
    config = Configuration()
    config.has_control = args.control is not None
    if args.linear_map is None:
        config.linear_map_name = args.graph_file_name.split(".")[0] + "_linear_map.npz"
    else:
        config.linear_map_name = args.linear_map
    config.fragment_length = int(args.fragment_length)
    config.read_length = int(args.read_length)
    if args.q_value_threshold is not None:
        config.q_value_threshold = float(args.q_value_threshold)
    return config


def get_intervals(args):
    iclass = UniqueIntervals  # Use Intervals to skip filter dup
    samples = iclass(parse_input_file(args.sample, args.graph))
    control_name = args.sample if args.control is None else args.control
    controls = iclass(parse_input_file(control_name, args.graph))


    return samples, controls


def get_callpeaks(args):
    config = get_confiugration(args)
    find_or_create_linear_map(args.graph, config.linear_map_name)
    out_name = args.out_name if args.out_name is not None else ""
    reporter = Reporter(out_name)
    return CallPeaks(args.graph, config, reporter)


def parse_input_file(input, graph):
    try:
        if isinstance(input, obg.IntervalCollection):
            return input
        elif input.endswith(".json"):
            intervals = vg_json_file_to_interval_collection(input, graph)
            return intervals
        else:
            try:
                intervals = obg.IntervalCollection.from_file(
                    input, graph=graph)
            except OSError:
                intervals = obg.IntervalCollection.from_file(
                    input, graph=graph, text_file=True)
            return intervals
    except FileNotFoundError as e:
        logging.debug(e)
        logging.critical("Input file %s not found. Aborting. " % input)
        sys.exit(1)

def run_callpeaks_interface(args):
    caller = get_callpeaks(args)
    inputs, controls = get_intervals(args)
    caller.run(inputs, controls)

    outname = args.out_name if args.out_name is not None else ""
    if args.sequence_graph is not None:
        PeakFasta(args.sequence_graph).write_max_path_sequences(
            outname+"sequences.fasta", caller.max_path_peaks)
    else:
        logging.info(
            "Not saving max path sequences, since a sequence graph was not found.")


def find_or_create_linear_map(graph, linear_map_name):

    if os.path.isfile(linear_map_name):
        logging.info("Found linear map %s. "
                     "Will be used in peak calling." % linear_map_name)
    else:
        logging.info("No linear map provided, and none found. Will create now.")
        create_linear_map(graph, linear_map_name)
        logging.info("Done creating linear map")


def run_callpeaks_whole_genome(args):
    logging.info("Running run_callpeaks_whole_genome")

    config = Configuration()
    if args.keep_duplicates == "True":
        config.keep_duplicates = True
        logging.info("Keeping duplicates")

    logging.info("Running whole genome.")
    chromosomes = args.chromosomes.split(",")
    graph_file_names = [args.data_dir + "/" + chrom + ".nobg" for chrom in chromosomes]
    logging.info("Will use graphs: %s" % graph_file_names)

    linear_map_file_names = []
    for i, chrom in enumerate(chromosomes):
        linear_map_name = args.data_dir + "/" + chrom + "_linear_map.npz"
        if not os.path.isfile(linear_map_name):
            logging.warning("Did not find linear map for "
                            "chromosome %s. Will create." % chrom)
            graph = obg.Graph.from_numpy_file(graph_file_names[i])
            create_linear_map(graph, linear_map_name)
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

    logging.info("Will use input alignments from %s" % sample_file_names)

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
        logging.info("Using input alignments as control")
        control_file_names = sample_file_names.copy()

    config.fragment_length = int(args.fragment_length)
    config.read_length = int(args.read_length)

    if config.fragment_length < config.read_length:
        logging.critical("Fragment length is smaller than read length. Cannot call peaks.")
        sys.exit(1)


    genome_size = int(args.genome_size)
    config.min_background = int(args.unique_reads) * int(args.fragment_length) / genome_size
    logging.info(
        "Computed min background signal to be %.3f using fragment length %f, "
        " %d unique reads, and genome size %d" % (config.min_background,
                                                  config.fragment_length,
                                                  int(args.unique_reads),
                                                  int(genome_size)))

    out_name = args.out_name if args.out_name is not None else ""
    reporter = Reporter(out_name)
    config.has_control = args.control is not None
    caller = MultipleGraphsCallpeaks(
        chromosomes,
        graph_file_names,
        sample_file_names,
        control_file_names,
        linear_map_file_names,
        config, reporter,
        sequence_retrievers=sequence_retrievers,
        stop_after_p_values=args.stop_after_p_values == "True",
    )
    caller.run()


def run_callpeaks_whole_genome_from_p_values(args):
    logging.info("Running whole genome from p-values.")
    chromosome = args.chromosome
    chromosomes = [chromosome]
    graph_file_names = [args.data_dir + chrom + ".nobg" for chrom in chromosomes]
    sequence_retrievers = \
        (obg.SequenceGraph.from_file(args.data_dir + "/" + chrom + ".nobg.sequences") for chrom in chromosomes)
    out_name = args.out_name if args.out_name is not None else ""
    config = Configuration()
    config.fragment_length = int(args.fragment_length)
    config.read_length = int(args.read_length)
    reporter = Reporter(out_name)
    caller = MultipleGraphsCallpeaks(
        chromosomes,
        graph_file_names,
        None,
        None,
        None,
        config,
        reporter,
        sequence_retrievers=sequence_retrievers,
    )
    caller.create_joined_q_value_mapping()
    caller.run_from_p_values(only_chromosome=chromosome)
