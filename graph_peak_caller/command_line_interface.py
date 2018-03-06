#!/usr/bin/python3
import argparse
import logging
import sys
import matplotlib as mpl
mpl.use('Agg')

import offsetbasedgraph as obg
from graph_peak_caller.peakcollection import Peak, PeakCollection

from graph_peak_caller.callpeaks_interface import \
    run_callpeaks_interface, run_callpeaks_whole_genome,\
    run_callpeaks_whole_genome_from_p_values

from graph_peak_caller.analysis_interface import analyse_peaks_whole_genome,\
    analyse_peaks, intervals_to_fasta, linear_peaks_to_fasta,\
    analyse_manually_classified_peaks, differential_expression,\
    plot_motif_enrichment

from graph_peak_caller.preprocess_interface import shift_estimation,\
    count_unique_reads_interface, create_ob_graph,\
    create_linear_map_interface,\
    split_vg_json_reads_into_chromosomes


logging.basicConfig(
    stream=sys.stdout, level=logging.WARNING,
    format="%(asctime)s, %(levelname)s: %(message)s")


def main():
    create_argument_parser()

def version(args):
    print("Graph Peak Caller v1.0.2")

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


def find_linear_path(args):
    from graph_peak_caller.util import create_linear_path
    import pyvg
    graph = obg.GraphWithReversals.from_numpy_file(args.ob_graph_file_name)
    vg_graph = pyvg.Graph.from_file(args.vg_json_graph_file_name)
    linear_path = create_linear_path(graph, vg_graph,
                                     path_name=args.linear_path_name,
                                     write_to_file=None)
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


def peaks_to_linear(args):
    # Get approximate linear position of peaks using a linear path
    linear_path = obg.NumpyIndexedInterval.from_file(args.linear_path_file_name)
    peaks = PeakCollection.from_file(args.peaks_file_name, text_file=True)
    linear_peaks = peaks.to_approx_linear_peaks(linear_path, args.chromosome)
    linear_peaks.to_bed_file(args.out_file_name)


interface = \
{
    'callpeaks':
        {
            'help': 'Callpeaks',
            'arguments':
                [
                    ('graph_file_name', ""),
                    ('vg_graph_file_name', "Graph file name (.vg). Used for fetching sequences. Can be set to 'None'."),
                    ('linear_map_base_name', "Set to base name of linear map."),
                    ('sample_reads_file_name', 'File name to a vg JSON file or intervalcollection file.'),
                    ('control_reads_file_name', 'File name to a vg JSON file or intervalcollection file. Set to the same as sample if no control.'),
                    ('with_control', 'True/False. Set to False if control is the same as sample.'),
                    ('out_base_name', 'eg experiment1_. Will be preprended to all output files.'),
                    ('fragment_length', 'The fragment length used in this ChIP-seq experiment. If unknown, set to an '
                                        'arbitrary number, e.g. 200. However, for good results, this number should be accurate.'),
                    ('read_length', 'The read length.')
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
                    ('sample_reads_base_name', 'Will use files *_[chromosome].json where * is the base name'),
                    ('control_reads_base_name', 'Will use files *_[chromosome].json where * is the base name'),
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
                    ('fimo_file_name', ''),
                    ('vg_graph_name', '')
                ],
            'method': differential_expression
        },
    'peaks_to_linear':
        {
            'help': 'Converts graph peaks to linear peaks by using a linear path through the graph',
            'arguments':
                [
                    ('peaks_file_name', "Shuold be a JSON intervalcollection, e.g. one created by callpeaks."),
                    ('linear_path_file_name', "Name of linera path file"),
                    ('chromosome', 'Name of chromosome that will be used when writing the bed file'),
                    ('out_file_name', 'Out file name')
                ],
            'method': peaks_to_linear
        },
    'version':
        {
            'help': 'Prints the current version',
            'arguments': [],
            'method': version
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

        subparser = subparsers.add_parser(
            command,
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
