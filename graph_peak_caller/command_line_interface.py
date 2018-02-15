#!/usr/bin/python3
import argparse
import logging
import os
import pickle
import sys

import offsetbasedgraph as obg
from pyvg.protoparser import json_file_to_obg_graph,\
    json_file_to_obg_numpy_graph
from pyvg.sequences import SequenceRetriever
from pyvg.util import vg_gam_file_to_interval_collection,\
    vg_json_file_to_interval_collection

from graph_peak_caller.callpeaks import CallPeaks, ExperimentInfo,\
    CallPeaksFromQvalues, Configuration
from graph_peak_caller.peakcollection import Peak
from graph_peak_caller.sparsepileup import SparsePileup
from graph_peak_caller.util import create_linear_map, create_ob_graph_from_vg
# from memory_profiler import profile
from graph_peak_caller.multiplegraphscallpeaks import MultipleGraphsCallpeaks
from graph_peak_caller.shift_estimation_multigraph import MultiGraphShiftEstimator
from pyvg.util import vg_gam_file_to_intervals

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s, %(levelname)s: %(message)s")

def main():
    create_argument_parser()

def shift_estimation(args):

    """
    from graph_peak_caller.shiftestimation import Treatment, Opt, PeakModel
    treatment = Treatment.from_bed_file("linear_bed.bed")
    opt = Opt()
    opt.gsize = 48172484
    model = PeakModel(opt, treatment)
    print(model.d)
    return
    """

    chromosomes = args.chromosomes.split(",")
    graphs = [args.ob_graphs_location + chrom for chrom in chromosomes]
    logging.info("Will try to use graphs %s" % graphs)
    sample_file_names = [args.sample_reads_base_name + chrom + ".json"
                         for chrom in chromosomes]
    logging.info("Will use reads from %s" % sample_file_names)

    estimator = MultiGraphShiftEstimator.from_files(
        chromosomes, graphs, sample_file_names)

    estimator.to_linear_bed_file("linear_bed.bed", read_length=36)

    d = estimator.get_estimates()
    print("Shift: %d" % d)


def filter_reads_in_graph(args):
    graph = obg.GraphWithReversals.from_unknown_file_format(args.graph_base_name)
    intervals = vg_gam_file_to_intervals(None, args.gam_reads_file, graph)

    filtered = []
    kept = 0
    i = 0
    print("Starting")
    for interval in intervals:
        if True or i % 10 == 0:
            logging.info("Processed %d intervals" % i)
        i += 1

        if interval.region_paths[0] in graph.blocks \
                and interval.region_paths[-1] in graph.blocks:
            filtered.append(interval)
            kept += 1
    logging.info("Kept %d intevals" % kept)

    filtered = obg.IntervalCollection(filtered)
    filtered.to_file(args.out_file_name)
    logging.info("Wrote resulting intervals to %s" % args.out_file_name)


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
             None, sample_file_name, ob_graph)
        control_intervals = vg_json_file_to_interval_collection(
             None, control_file_name, ob_graph)
    else:
        sample_intervals = obg.IntervalCollection.from_file(
            sample_file_name, graph=ob_graph)
        control_intervals = obg.IntervalCollection.from_file(
            control_file_name, graph=ob_graph)

    graph_size = ob_graph.number_of_basepairs()
    logging.info("Number of base pairs in graph: %d" % graph_size)

    experiment_info = ExperimentInfo(graph_size, fragment_length, read_length)
    config = Configuration(
        skip_read_validation=True, save_tmp_results_to_file=False,
        skip_filter_duplicates=False, p_val_cutoff=0.05,
        graph_is_partially_ordered=True)

    caller = CallPeaks.run_from_intervals(
        ob_graph, sample_intervals, control_intervals,
        experiment_info=experiment_info,
        out_file_base_name=out_name, has_control=has_control,
        linear_map=linear_map_file_name,
        configuration=config
    )
    retriever = SequenceRetriever.from_vg_graph(vg_graph_file_name)
    caller.save_max_path_sequences_to_fasta_file("sequences.fasta", retriever)


def intervals_to_fasta(args):
    logging.info("Getting sequence retriever")
    retriever = SequenceRetriever.from_vg_graph(args.vg_graph_file_name)
    logging.info("Getting intervals")
    intervals = obg.IntervalCollection.create_generator_from_file(
        args.intervals_file_name)
    logging.info("Writing to fasta")
    CallPeaksFromQvalues.intervals_to_fasta_file(
        intervals, args.out_file_name, retriever)


# @profile
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

def count_unique_reads_interface(args):
    chromosomes = args.chromosomes.split(",")
    graph_file_names = [args.graphs_location + chrom for chrom in chromosomes]
    reads_file_names = [args.reads_base_name + chrom + ".json"
                        for chrom in chromosomes]

    count_unique_reads(chromosomes, graph_file_names,
                                  reads_file_names)


def count_unique_reads(chromosomes, graph_file_names, reads_file_names):

    graphs = (obg.GraphWithReversals.from_numpy_file(f) for f in graph_file_names)
    reads = (vg_json_file_to_interval_collection(None, f, graph)
             for f, graph in zip(reads_file_names, graphs))

    unique_reads = MultipleGraphsCallpeaks.count_number_of_unique_reads(reads)
    print(unique_reads)


def run_callpeaks_whole_genome(args):
    logging.info("Running whole genome.")
    chromosomes = args.chromosomes.split(",")
    graph_file_names = [args.graphs_location + chrom for chrom in chromosomes]
    linear_map_file_names = [args.linear_maps_location + chrom for chrom in chromosomes]
    vg_graphs = [args.vg_graphs_location + chrom + ".vg" for chrom in chromosomes]
    sequence_retrievers = \
        (SequenceRetriever.from_vg_graph(fn) for fn in vg_graphs)
    sample_file_names = [args.sample_reads_base_name + chrom + ".json"
                        for chrom in chromosomes]
    control_file_names = [args.sample_reads_base_name + chrom + ".json"
                        for chrom in chromosomes]


    min_background = int(args.unique_reads) * int(args.fragment_length) / int(args.genome_size)
    logging.info("Computed min background signal to be %.3f" % min_background)

    caller = MultipleGraphsCallpeaks(
        chromosomes,
        graph_file_names,
        sample_file_names,
        control_file_names,
        linear_map_file_names,
        int(args.fragment_length),
        int(args.read_length),
        has_control=args.with_control=="True",
        sequence_retrievers=sequence_retrievers,
        out_base_name=args.out_base_name,
        stop_after_p_values=args.stop_after_p_values == "True",
        min_background=min_background
    )
    caller.run()

def run_callpeaks_whole_genome_from_p_values(args):
    logging.info("Running whole genome from p-values.")
    chromosome = args.chromosome
    chromosomes = args.chromosomes.split(",")
    graph_file_names = [args.graphs_location + chrom for chrom in chromosomes]
    linear_map_file_names = [args.linear_maps_location + chrom for chrom in chromosomes]
    vg_graphs = [args.vg_graphs_location + chrom + ".vg" for chrom in chromosomes]
    sequence_retrievers = \
        (SequenceRetriever.from_vg_graph(fn) for fn in vg_graphs)
    sample_file_names = [args.sample_reads_base_name + chrom + ".json"
                        for chrom in chromosomes]
    control_file_names = [args.sample_reads_base_name + chrom + ".json"
                        for chrom in chromosomes]

    caller = MultipleGraphsCallpeaks(
        chromosomes,
        graph_file_names,
        sample_file_names,
        control_file_names,
        linear_map_file_names,
        int(args.fragment_length),
        int(args.read_length),
        has_control=args.with_control=="True",
        sequence_retrievers=sequence_retrievers,
        out_base_name=args.out_base_name
    )
    caller.create_joined_q_value_mapping()
    caller.run_from_p_values(only_chromosome=chromosome)

def linear_peaks_to_fasta(args):
    from graph_peak_caller.nongraphpeaks import NonGraphPeakCollection
    collection = NonGraphPeakCollection.from_bed_file(args.linear_reads_file_name)
    collection.set_peak_sequences_using_fasta(fasta_file_location=args.fasta_file)
    collection.save_to_sorted_fasta(args.out_file_name)
    logging.info("Saved sequences to %s" % args.out_file_name)

def create_ob_graph(args):
    logging.info("Creating obgraph")
    ob_graph = json_file_to_obg_numpy_graph(args.vg_json_file_name, 0)

    logging.info("Writing ob graph to file")
    ob_graph.to_numpy_file(args.out_file_name)


def create_linear_map_interface(args):
    logging.info("Reading ob graph from file")
    ob_graph = obg.GraphWithReversals.from_numpy_file(args.obg_file_name)
    logging.info("Converting to dict format, allowing graph to "
                 "be changed by linear map process")
    ob_graph.convert_to_dict_backend()
    logging.info("Creating linear map")
    create_linear_map(ob_graph, args.vg_snarls_file_name, args.out_file_base_name, copy_graph=False)


def split_vg_json_reads_into_chromosomes(args):
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
        print("Processing chromosome %s" % chromosome)
        fasta_file = open(chromosome + "_sequences.fasta")
        for line in fasta_file:
            if line.startswith(">"):
                all_fasta_entries.append([line, None])
            else:
                # This is sequence, add to prev entry
                all_fasta_entries[-1][1] = line
        #fasta_file.close()

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
    from benchmarking.motifenrichment import plot_true_positives
    fasta1 = args.fasta1
    fasta2 = args.fasta2
    meme = args.meme_motif_file

    plot_true_positives(
        {
            "file1": fasta1,
            "file2": fasta2
        },
        meme,
        save_to_file=args.out_figure_file_name,
        run_fimo=args.run_fimo == "True"
    )

def analyse_peaks(args):
    from graph_peak_caller.peakscomparer import PeaksComparerV2
    from graph_peak_caller.analyse_peaks import LinearRegion
    graph = obg.GraphWithReversals.from_numpy_file(args.ob_graph_file_name)

    end = int(args.graph_end)
    if end == 0:
        region = None
    else:
        region = LinearRegion(args.graph_chromosome, int(args.graph_start), end)
    analyser = PeaksComparerV2(graph,
                               args.vg_graph_file_name,
                               args.linear_peaks_fasta_file_name,
                               args.graph_peaks_fasta_file_name,
                               args.linear_peaks_fimo_results_file,
                               args.graph_peaks_fimo_results_file,
                               linear_path_name=args.linear_path_name,
                               region=region)


def analyse_peaks_whole_genome(args):
    from graph_peak_caller.peakscomparer import PeaksComparerV2, AnalysisResults
    from offsetbasedgraph import IndexedInterval
    chromosomes = args.chromosomes.split(",")
    graph_file_names = (args.graphs_dir + chrom + ".nobg" for chrom in chromosomes)
    graphs = (obg.GraphWithReversals.from_numpy_file(fn) for fn in graph_file_names)
    vg_file_names = (args.graphs_dir + chrom + ".vg" for chrom in chromosomes)


    results = AnalysisResults()

    for chrom in chromosomes:
        graph = obg.GraphWithReversals.from_numpy_file(
            args.graphs_dir + chrom + ".nobg")
        logging.info("Reading linear path")
        linear_path = IndexedInterval.from_file(args.graphs_dir + chrom + "_linear_path.interval")

        analyser = PeaksComparerV2(
            graph,
            args.graphs_dir + chrom + ".json",
            args.results_dir + "macs_sequences_chr%s.fasta" % chrom,
            args.results_dir + "%s_sequences.fasta" % chrom,
            args.results_dir + "/fimo_macs_chr%s/fimo.txt" % chrom,
            args.results_dir + "/fimo_graph_chr%s/fimo.txt" % chrom,
            linear_path
        )
        results = results + analyser.results

    print(" === Final results for all chromosomes ===")
    print(results)

    results.to_file(args.out_file + ".pickle")
    with open(args.out_file + ".txt", "w") as f:
        f.write(str(results))
    logging.info("Wrote results as pickle to %s and as text to %s" \
                 % (args.out_file + ".pickle", args.out_file + ".txt"))


def analyse_manually_classified_peaks(args):
    for chromosome in args.chromosomes.split():
        graph_file_name = args.graphs_location + "/" + chromosome + ".nobg"
        graph = obg.GraphWithReversals.from_numpy_file(graph_file_name)

        from .manually_classified_peaks import CheckOverlapWithManuallyClassifiedPeaks
        CheckOverlapWithManuallyClassifiedPeaks.from_graph_peaks_in_fasta(
                graph,
                args.graphs_location + "/" + chromosome  + ".json",
                chromosome,
                args.reads_base_name + chromosome + "_sequences.fasta",
                args.regions_file,
                args.manually_classified_peaks_file)


def find_linear_path(args):
    from graph_peak_caller.util import create_linear_path
    import pyvg
    graph = obg.GraphWithReversals.from_numpy_file(args.ob_graph_file_name)
    vg_graph = pyvg.vg.Graph.create_from_file(args.vg_json_graph_file_name)
    linear_path = create_linear_path(graph, vg_graph, path_name=args.linear_path_name)
    linear_path.to_file(args.out_file_name)
    logging.info("Wrote to file %s" % args.out_file_name)


interface = \
{

    'callpeaks':
        {
            'help': 'Callpeaks',
            'arguments':
                [
                    ('graph_file_name', ""),
                    ('vg_graph_file_name', "Graph file name (.vg)"),
                    ('linear_map_base_name', "Set to desired base name. Will be used if exists, created if not."),
                    ('sample_reads_file_name', ' '),
                    ('control_reads_file_name', ' '),
                    ('with_control', 'True/False'),
                    ('out_base_name', 'eg experiment1_'),
                    ('fragment_length', ''),
                    ('read_length', '')
                ],
            'method': run_callpeaks_interface
        },
    'callpeaks_whole_genome':
        {
            'help': 'Callpeaks on whole genome, using one graph for each chromosome',
            'arguments':
                [
                    ('chromosomes', 'Comma-separated list of chromosomes to use, e.g. 1,2,X,8,Y'),
                    ('graphs_location', 'Will use the graphs *_[chromosome]'),
                    ('vg_graphs_location', ''),
                    ('linear_maps_location', ''),
                    ('sample_reads_base_name', 'Will use files *_[chromosome].json'),
                    ('control_reads_base_name', 'Will use files *_[chromosome].json'),
                    ('out_base_name', 'eg experiment1_'),
                    ('with_control', 'True/False'),
                    ('fragment_length', ''),
                    ('read_length', ''),
                    ('stop_after_p_values', 'True/False - whether to only run until p-value track is computed (before peak calling)'),
                    ('unique_reads', 'Number of unique reads. Found by calling count_unique_reads'),
                    ('genome_size', 'Number of base pairs covered by graphs in total (on a linear genome)')
                ],
            'method': run_callpeaks_whole_genome
        },
    'callpeaks_whole_genome_from_p_values':
        {
            'help': 'Callpeaks on whole genome, using one graph for each chromosome',
            'arguments':
                [
                    ('chromosome', 'Specific chromosome'),
                    ('chromosomes', 'All chromosomes to compute q-values from, e.g. 1,2,3,4,X,Y'),
                    ('graphs_location', 'Will use the graphs *_[chromosome]'),
                    ('vg_graphs_location', ''),
                    ('linear_maps_location', ''),
                    ('sample_reads_base_name', 'Will use files *_[chromosome].json'),
                    ('control_reads_base_name', 'Will use files *_[chromosome].json'),
                    ('out_base_name', 'eg experiment1_'),
                    ('with_control', 'True/False'),
                    ('fragment_length', ''),
                    ('read_length', '')
                ],
            'method': run_callpeaks_whole_genome_from_p_values
        },
    'estimate_shift':
        {
            'help': 'Estimate shift on whole genome',
            'arguments':
                [
                    ('chromosomes', 'Comma-separated list of chromosomes to use, e.g. 1,2,X,8,Y'),
                    ('ob_graphs_location', 'Location of graph files'),
                    ('sample_reads_base_name', 'Will use files [base_name][chromosome].json'),
                ],
            'method': shift_estimation
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
            'help': 'Creates Offset Based Graph from a vg json graph file. Stores the resulting graph to a file.',
            'arguments':
                [
                    ('vg_json_file_name', 'Vg json file name (created by running vg view -Vj graph.vg > graph.json'),
                    ('out_file_name', 'E.g. graph.obg')
                ],
            'method': create_ob_graph
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
                    ('out_figure_file_name', ''),
                    ('run_fimo', 'Set to True if fimo has not already been run.')
                ],
            'method': plot_motif_enrichment
        },
    'filter_reads_in_graph':
        {
            'help': "Writes reads in graph to new json file.",
            'arguments':
                [
                    ('gam_reads_file', ''),
                    ('graph_base_name', ''),
                    ('out_file_name', '')
                ],
            'method': filter_reads_in_graph
        },
    'analyse_peaks':
        {
            'help': 'Analyse linear peaks and graph peaks.',
            'arguments':
                [
                    ('ob_graph_file_name', ''),
                    ('vg_graph_file_name', ''),
                    ('linear_peaks_fasta_file_name', ''),
                    ('graph_peaks_fasta_file_name', ''),
                    ('linear_peaks_fimo_results_file', ''),
                    ('graph_peaks_fimo_results_file', ''),
                    ('linear_path_name', 'Ref of chromosome name. Used to find path in vg json file'),
                    ('graph_chromosome', 'None if not relevant'),
                    ('graph_start', 'Start pos in chromosome of graph.'),
                    ('graph_end', 'End pos in chromosome of graph. 0 if covering whole chromosome')
                ],
            'method': analyse_peaks
        },
    'analyse_peaks_whole_genome':
        {
            'help': "Analyse peak results from whole genome run",
            'arguments':
                [
                    ('chromosomes', 'Comma seaparated list of chromosomes to use'),
                    ('results_dir', 'Directory of result files (should contain fasta files from both linear and graph peak calling'),
                    ('graphs_dir', 'Dir containing obg graphs on form 1.nobg, 2.nobg, ... and vg graphs 1.json, 2.json, ...'),
                    ('out_file', 'Out file base name (file endings for different formats will be appended)')
                ],
            'method': analyse_peaks_whole_genome
        },
    'count_unique_reads':
        {
            'help': 'Count unique reads for whole genome set of read files.',
            'arguments':
                [
                    ('chromosomes', 'Comma-separated list of chromosomes to use, e.g. 1,2,X,8,Y'),
                    ('graphs_location', 'Will use the graphs *_[chromosome]'),
                    ('reads_base_name', 'Will use files *_[chromosome].json')
                ],
            'method': count_unique_reads_interface
        },
    'analyse_manually_classified_peaks':
        {
            'help': '',
            'arguments':
                [
                    ('chromosomes', 'Comma-separated list of chromosomes to use, e.g. 1,2,X,8,Y'),
                    ('graphs_location', 'Will use the graphs *_[chromosome]'),
                    ('reads_base_name', 'Will use files *_[chromosome].json'),
                    ('regions_file', ''),
                    ('manually_classified_peaks_file', '')
                ],
            'method': analyse_manually_classified_peaks
        },
    'find_linear_path':
        {
            'help': 'Finds lineat path through graph. Saves as indexed interval to file.',
            'arguments':
                [
                    ('vg_json_graph_file_name', ''),
                    ('ob_graph_file_name', ''),
                    ('linear_path_name', 'Name of path in the vg graph (typically ref or chromosome name'),
                    ('out_file_name', ''),
                ],
            'method': find_linear_path
        }
}




def create_argument_parser():
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

if __name__ == "__main__":
    create_argument_parser()
