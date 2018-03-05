#!/usr/bin/python3
import argparse
import logging
import sys

import offsetbasedgraph as obg
from pyvg.sequences import SequenceRetriever
from pyvg.conversion import vg_json_file_to_interval_collection,\
    json_file_to_obg_numpy_graph

from graph_peak_caller import CallPeaksFromQvalues
from graph_peak_caller.peakcollection import Peak
from graph_peak_caller.util import create_linear_map
from graph_peak_caller.multiplegraphscallpeaks import MultipleGraphsCallpeaks
from graph_peak_caller.shift_estimation_multigraph import \
    MultiGraphShiftEstimator

from graph_peak_caller.callpeaks_interface import \
    run_callpeaks_interface, run_callpeaks_whole_genome,\
    run_callpeaks_whole_genome_from_p_values

from graph_peak_caller.analysis_interface import analyse_peaks_whole_genome,\
    analyse_peaks,\
    analyse_manually_classified_peaks, differential_expression

logging.basicConfig(
    stream=sys.stdout, level=logging.WARNING,
    format="%(asctime)s, %(levelname)s: %(message)s")


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


def intervals_to_fasta(args):
    logging.info("Getting sequence retriever")
    retriever = SequenceRetriever.from_vg_graph(args.vg_graph_file_name)
    logging.info("Getting intervals")
    intervals = obg.IntervalCollection.create_generator_from_file(
        args.intervals_file_name)
    logging.info("Writing to fasta")
    CallPeaksFromQvalues.intervals_to_fasta_file(
        intervals, args.out_file_name, retriever)


def count_unique_reads_interface(args):
    chromosomes = args.chromosomes.split(",")
    graph_file_names = [args.graphs_location + chrom for chrom in chromosomes]
    reads_file_names = [args.reads_base_name + chrom + ".json"
                        for chrom in chromosomes]

    count_unique_reads(chromosomes, graph_file_names,
                       reads_file_names)


def count_unique_reads(chromosomes, graph_file_names, reads_file_names):

    graphs = (obg.GraphWithReversals.from_numpy_file(f)
              for f in graph_file_names)
    reads = (vg_json_file_to_interval_collection(f, graph)
             for f, graph in zip(reads_file_names, graphs))

    unique_reads = MultipleGraphsCallpeaks.count_number_of_unique_reads(reads)
    print(unique_reads)


def linear_peaks_to_fasta(args):
    from graph_peak_caller.nongraphpeaks import NonGraphPeakCollection
    collection = NonGraphPeakCollection.from_bed_file(
        args.linear_reads_file_name)
    collection.set_peak_sequences_using_fasta(
        fasta_file_location=args.fasta_file)
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
    create_linear_map(ob_graph, args.vg_snarls_file_name,
                      args.out_file_base_name, copy_graph=False)


def split_vg_json_reads_into_chromosomes(args):
    reads_base_name = args.vg_json_reads_file_name.split(".")[0]
    logging.info("Will write reads to files %s_[chromosome].json" % reads_base_name)

    chromosomes = args.chromosomes.split(",")
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
        mapped_chrom = None
        for chrom in chromosomes:
            if node >= chromosome_limits[chrom][0] and node <= chromosome_limits[chrom][1]:
                mapped_chrom = chrom
                break
        #assert mapped_chrom is not None, "Found no match for node id %d" % node
        return mapped_chrom

    n_without_node_id = 0
    for line in reads_file:
        if i % 100000 == 0:
            logging.info("Line #%d" % i)
        i += 1

        groups = regex.search(line)
        if groups is None:
            n_without_node_id += 1
            continue
        groups = groups.groups()
        if len(groups) > 0:
            node = int(groups[0])
            mapped_chrom = get_mapped_chrom(node)
            if mapped_chrom is None:
                n_without_node_id += 1
                continue
            out_files[mapped_chrom].writelines([line])
        else:
            print("No groups fond")

    for file in out_files.values():
        file.close()

    logging.info("Done. Found %d lines without node id or not matching into the given list of chromosomes" % n_without_node_id)


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
    from graph_peak_caller.motifenrichment import plot_true_positives
    fasta1 = args.fasta1
    fasta2 = args.fasta2
    meme = args.meme_motif_file

    plot_true_positives(
        [
            ("Graph Peak Caller", fasta1),
            ("MACS2", fasta2)
        ],
        meme,
        plot_title=args.plot_title,
        save_to_file=args.out_figure_file_name,
        run_fimo=args.run_fimo == "True"
    )


def find_linear_path(args):
    from graph_peak_caller.util import create_linear_path
    import pyvg
    graph = obg.GraphWithReversals.from_numpy_file(args.ob_graph_file_name)
    vg_graph = pyvg.Graph.from_file(args.vg_json_graph_file_name)
    linear_path = create_linear_path(graph, vg_graph,
                                     path_name=args.linear_path_name, write_to_file=None)
    linear_path.to_file(args.out_file_name)
    logging.info("Wrote to file %s" % args.out_file_name)


def move_linear_reads_to_graph(args):
    chromosomes = args.chromosomes.split(",")
    chrom_lookup = set(chromosomes)
    graphs = {}
    out_files = {}
    linear_paths = {}
    for chrom in chromosomes:
        linear_paths[chrom] = obg.NumpyIndexedInterval.from_file(
            args.data_dir + "/" + chrom + "_linear_pathv2.interval"
        )
        graphs[chrom] = obg.Graph.from_numpy_file(args.data_dir + "/" + chrom + ".nobg")
        out_files[chrom] = open(args.out_files_base_name + "_" + chrom + ".intervalcollection", "w")

    bedfile = open(args.bed_file_name, "r")
    i = 0
    for line in bedfile:
        if i % 100000 == 0:
            logging.info("%d reads processed" % i)
        i += 1
        line = line.split()
        is_reverse = line[5] == "-"
        chrom = line[0].replace("chr", "")
        if chrom not in chrom_lookup:
            continue
        start = int(line[1])
        end = int(line[2])
        graph_interval = linear_paths[chrom].get_exact_subinterval(start, end)
        graph_interval.graph = graphs[chrom]

        if is_reverse:
            graph_interval = graph_interval.get_reverse()
            assert graph_interval.region_paths[0] < 0
        out_files[chrom].writelines(["%s\n" % graph_interval.to_file_line()])

    for chrom in chromosomes:
        out_files[chrom].close()


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
            'help': 'Callpeaks on whole genome from pvalues. Assumes pvalue files are already computed.',
            'arguments':
                [
                    ('chromosome', 'Specific chromosome to find peaks for.'),
                    ('graphs_location', 'Directory containing graphs.'),
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
                    ('chromosomes', 'Comma-separated list of chromosomes to split reads into'),
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
                    ('run_fimo', 'Set to True if fimo has not already been run.'),
                    ('plot_title', 'Title above plot')
                ],
            'method': plot_motif_enrichment
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
        },
    'move_linear_reads_to_graph':
        {
            'help': 'Translate reads in sam file to reads on a graph. Writes to interval collection files.',
            'arguments':
                [
                    ('bed_file_name', ''),
                    ('chromosomes', 'Comma separated list of chromosomes to use'),
                    ('data_dir', 'Directory containing graphs and linear path files'),
                    ('out_files_base_name', '')
                ],
            'method': move_linear_reads_to_graph
        },
    'diffexpr':
        {
            'help': 'Find differentially expressed motif matches',
            'arguments':
                [
                    ('test_name', ''),
                    ('graph_name', ''),
                    ('vg_graph_name', '')
                ],
            'method': differential_expression
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
