#!/usr/bin/python3
import argparse
import sys
import logging
from graph_peak_caller.logging_config import set_logging_config
from graph_peak_caller.custom_exceptions import *
import pyvg

import matplotlib as mpl
mpl.use('Agg')  # Required for server usage (e.g. travis)

import offsetbasedgraph as obg
from graph_peak_caller.peakcollection import Peak, PeakCollection
from graph_peak_caller.sparsediffs import SparseValues
from graph_peak_caller.mindense import DensePileup

from graph_peak_caller.callpeaks_interface import \
    run_callpeaks_interface, run_callpeaks_whole_genome,\
    run_callpeaks_whole_genome_from_p_values, run_callpeaks2

from graph_peak_caller.analysis.analysis_interface import analyse_peaks_whole_genome,\
    analyse_peaks, peaks_to_fasta, linear_peaks_to_fasta,\
    analyse_manually_classified_peaks, differential_expression,\
    plot_motif_enrichment, get_summits, peaks_to_linear, move_linear_reads_to_graph,\
    find_linear_path, concatenate_sequence_files, check_haplotype, get_motif_locations

from graph_peak_caller.preprocess_interface import \
    count_unique_reads_interface, create_ob_graph,\
    create_linear_map_interface,\
    split_vg_json_reads_into_chromosomes, shift_estimation


from offsetbasedgraph.interval import NoLinearProjectionException

def main():
    run_argument_parser(sys.argv[1:])


def version(args):
    print("Graph Peak Caller v1.1.1")


def project_alignments(alignments, linear_path):
    import numpy as np
    i = 0
    for alignment in alignments:
        if i % 10000 == 0:
            logging.info("Processed %d alignments" % i)
        i += 1

        try:
            if alignment.region_paths[0] < 0:
                assert np.all(np.array(alignment.region_paths) < 0)
                start, end = alignment.get_reverse().to_linear_offsets2(linear_path)
                strand = "-"
            else:
                assert np.all(np.array(alignment.region_paths) > 0)
                start, end = alignment.to_linear_offsets2(linear_path)
                strand = "+"
        except NoLinearProjectionException:
            logging.warning("Found no linear projection for %s. Skipping" % alignment)

        yield start, end, strand

def project_vg_alignments(args):
    from pyvg.conversion import vg_json_file_to_intervals
    linear_path = obg.NumpyIndexedInterval.from_file(args.linear_path_file_name)
    #print(linear_path.nodes_in_interval())
    alignments = vg_json_file_to_intervals(args.alignments_json_file_name, args.graph)

    with open(args.out_file_name, "w") as f:
        projected_alignments = project_alignments(alignments, linear_path)
        for start, end, strand in projected_alignments:
            f.writelines(["%s\t%d\t%d\t.\t0\t%s\n" % (args.chromosome, start, end, strand)])


def check_pruned_graphs_stats(args):
    from pyvg.util import get_stats
    chromosomes = args.chromosomes.split(",")
    graphs_dir = args.graphs_dir

    nodes_pruned = {}
    lengths_pruned = {}
    edges_pruned = {}

    for chromosome in chromosomes:
        stats_pruned = get_stats(graphs_dir + "/" + chromosome + ".prune.vg")
        stats_nonpruned = get_stats(graphs_dir + "/" + chromosome + ".vg")

        print(stats_pruned)
        print(stats_nonpruned)


def vg_json_alignments_to_intervals(args):
    from pyvg.conversion import vg_json_file_to_interval_collection
    interval_collection = vg_json_file_to_interval_collection(args.vg_json_file_name, args.graph)
    interval_collection.to_file(args.out_file_name, text_file=True)
    logging.info("Wrote to file %s" % args.out_file_name)


def get_intersecting_intervals(args):
    from offsetbasedgraph import IntervalCollection
    intervals1 = IntervalCollection.from_file(args.file1, text_file=True, graph=args.graph)
    intervals2 = IntervalCollection.from_file(args.file2, text_file=True, graph=args.graph)

    out = []
    for interval1 in intervals1.intervals:
        for interval2 in intervals2.intervals:
            if interval1.intersects(interval2):
                out.append(interval1)
                logging.info("Found match between %s and %s" % (interval1, interval2))
                continue

    IntervalCollection(out).to_file(args.out_file_name, text_file=True)
    logging.info("Wrote intersecting intervals to %s" % args.out_file_name)



interface = \
{

    'callpeaks':
        {
            'help': 'Call peaks one or multiple graphs.',
            'arguments':
                [
                    ('-g/--graph', 'Graph(s) file name. Either a single file name, e.g. graph.nobg or a '
                                   'pattern with the * wildcard, such as graph_*.nobg. In the latter case'
                                   ', * will be replaced by the names specified in --chromosomes.'),
                    ('-s/--sample', 'Sample alignments. Use wildcard * to specify multiple files'),
                    ('-c/--control', 'Optional. Use wildcard * to specify multiple files.'
                                     ' Use only this option if you have separate control alignments'),
                    ('-f/--fragment_length', 'Optional. Specify if you know the exact fragment length. '
                                             'Will be estimated if not specified.'),
                    ('-r/--read_length', 'Optional. '),
                    ('-p/--stop_after_p_values', 'Optional. Default False. Set to True in order to'
                                                 'stop the whole run after p-values has been computed.'
                                                 ' Useful if wanting to run multiple chromosomes in parallell, and'
                                                 'waiting before continuing from p-values.'),
                    ('-n/--out_name', 'Optional. Out base name. Prepended to output files.'),
                    ('-u/--unique_reads', 'Optional. Number of unique reads. '
                                          'Found by calling count_unique_reads'),
                    ('-G/--genome_size', 'Optional. Number of base pairs covered by '
                                         'graphs in total (on a linear genome). '
                                         'If not set, will be estimated if running on single graph. '
                                         'Must be set if running on multiple graphs.'),
                    ('-D/--keep_duplicates', 'Optional. Set to True in order to keep '
                                             'duplicate input alignments.'),
                    ('-m/--min_fold_enrichment', 'Optional. Minimum fold enrichment required for '
                                               'candidate peaks when estimating fragment length. Default 5.'),
                    ('-M /--max_fold_enrichment', 'Optional. Maximum fold enrichment required for '
                                               'candidate peaks when estimating fragment length. Default 50.'),

                ],
                'method': run_callpeaks2,
        },

    'callpeaks_whole_genome_from_p_values':
        {
            'help': 'Callpeaks on whole genome from pvalues.',
            'arguments':
                [
                    ('chromosome', 'Specific chromosome to find peaks for.'),
                    ('-d/--data_dir', 'Directory containing graphs.'),
                    ('-n/--out_name', 'Optional. eg experiment1_'),
                    ('-f/--fragment_length', ''),
                    ('-r/--read_length', ''),
                    ('-q/--q_threshold', 'Optional. q-value threshold. Default is 0.05.'),
                ],
            'method': run_callpeaks_whole_genome_from_p_values
        },
    'peaks_to_fasta':
        {
            'help': 'Get sequences for intervals in interval file.',
            'arguments':
                [
                    ('sequence_graph', 'Sequence Graph, e.g. graph.nobg.sequences '
                                       '(created by calling create_ob_graph'),
                    ('intervals_file_name', ''),
                    ('out_file_name', '')
                ],
            'method': peaks_to_fasta
        },
    'linear_peaks_to_fasta':
        {
            'help': 'Converts bed file of peaks to fasta',
            'arguments':
                [
                    ('linear_reads_file_name', 'E.g a macs file, peaks.narrowPeak'),
                    ('fasta_file', 'Reference genome fasta file. Will be used to fetch sequences.'),
                    ('out_file_name', '')
                ],
            'method': linear_peaks_to_fasta
        },
    'linear_peaks_to_fasta_summits':
        {
            'help': 'Converts bed file of peaks to fasta (but keep only summits)',
            'arguments':
                [
                    ('linear_reads_file_name', 'E.g a macs file, peaks.narrowPeak'),
                    ('fasta_file', 'Reference genome fasta file. Will be used to fetch sequences.'),
                    ('out_file_name', ''),
                    ('window', 'Number of bps around summits to keep')
                ],
            'method': linear_peaks_to_fasta
        },
    'create_ob_graph':
        {
            'help': 'Create OffsetBased Graph from vg.',
            'arguments':
                [
                    ('vg_json_file_name', 'Vg json file name '
                                          '(created by running vg view -Vj graph.vg > graph.json'),
                    ('-o/--out_file_name', 'Optional. Will use input file base name if unset.')
                ],
            'method': create_ob_graph
        },
    'create_linear_map':
        {
            'help': 'Create linear map from vg snarls.',
            'requires_graph': True,
            'arguments':
                [
                    ('-o/--out_file_base_name', 'Optional out file name'),
                ],
            'method': create_linear_map_interface
        },
    'split_vg_json_reads_into_chromosomes':
        {
            'help': "Split vg json reads by chromosome.",
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
            'help': "Merge multiple *_sequence.fasta files.",
            'arguments':
                [
                    ('chromosomes', 'comma delimted, e.g 1,2,3, used to fetch files of type chr1_sequences.fasta, ...'),
                    ('out_file_name', ''),
                    ('-s/--is_summits', 'Optional. Set to True if input files are *_sequences_summits.fasta.'
                                        'If False or not set, will search for *_sequences.fasta'),
                    ('-f/--use_input_file_pattern', 'Optional. If set, will look for input files'
                                                    'following the pattern. Use [chrom] where chromosome should'
                                                    'be replaced. E.g. peaks_[chrom].fasta')
                ],
            'method': concatenate_sequence_files
        },
    'plot_motif_enrichment':
        {
            'help': "Plots motif enrichments using fimo.",
            'arguments':
                [
                    ('fasta1', ''),
                    ('fasta2', ''),
                    ('meme_motif_file', 'something.meme'),
                    ('out_figure_file_name', ''),
                    ('run_fimo', 'Set to True if fimo has not already been run.'),
                    ('plot_title', 'Title above plot'),
                    ('background_model_file', 'Optional. Background model file for fimo. If not set, fimo will use uniform background.')
                ],
            'method': plot_motif_enrichment
        },
    'analyse_peaks':
        {
            'help': 'Analyse linear peaks and graph peaks.',
            'requires_graph': True,
            'arguments':
                [
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
            'help': 'Count unique reads in vg json alignments.',
            'arguments':
                [
                    ('chromosomes', 'Comma-separated list of chromosomes to use, e.g. 1,2,X,8,Y'),
                    ('graphs_location', 'Will use the graphs *_[chromosome].nobg'),
                    ('reads_base_name', 'Will use files *_[chromosome].json')
                ],
            'method': count_unique_reads_interface
        },
    'find_linear_path':
        {
            'help': 'Find linear path through graph.',
            'requires_graph': True,
            'arguments':
                [
                    ('vg_json_graph_file_name', ''),
                    ('linear_path_name', 'Name of path in the vg graph (typically ref or chromosome name'),
                    ('out_file_name', ''),
                ],
            'method': find_linear_path
        },
    'move_linear_reads_to_graph':
        {
            'help': 'Convert SAM to Intervals on graph.',
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
            'requires_graph': True,
            'arguments':
                [
                    ('test_name', ''),
                    ('fimo_file_name', '')
                ],
            'method': differential_expression
        },
    'check_haplotype':
        {
            'help': 'Check motifs for haplotype',
            'arguments':
                [
                    ('data_folder', ''),
                    ('fasta_file', ''),
                    ('result_folder', ''),
                    ('chrom', ''),
                    ('interval_name', ''),
                    ("all_reads", "Optional, wheter to check all reads covering region. Default: True"),
                    ("strict", "Optional, wheter to consider '.'s in the vcf as reference. Default: True")
                ],
            'method': check_haplotype
        },
    'motif_locations':
        {
            'help': 'Check motifs for haplotype',
            'arguments':
                [
                    ('data_folder', ''),
                    ('result_folder', ''),
                    ('chrom', '')
                ],
            'method': get_motif_locations
        },

    'peaks_to_linear':
        {
            'help': 'Converts graph peaks to linear peaks.',
            'arguments':
                [
                    ('peaks_file_name', "Shuold be a JSON intervalcollection, e.g. one created by callpeaks."),
                    ('linear_path_file_name', "Name of linear path file"),
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
        },
    'estimate_shift':
        {
            'help': 'Estimate shift using one or multiple graphs.',
            'arguments':
                [
                    ('chromosomes', 'Graph base names. Set to empty string if only single '
                                    'graph is being used. If whole-genome, use comma-separated '
                                    'list of chromosomes to use, e.g. 1,2,X,8,Y'),
                    ('ob_graphs_location', 'Location of graph files'),
                    ('sample_reads_base_name', 'Will use files [base_name][chromosome].json'),
                    ('min_fold_enrichment', 'Optional. Minimum fold enrichment when '
                                            'finding candidates. Default is 5.'),
                    ('max_fold_enrichment', 'Optional. Maximum fold enrichment when '
                                            'finding candidates. Default is 50.')
                ],
            'method': shift_estimation
        },
    'get_summits':
        {
            'help': 'Get summit around peaks. Will write to a new file, using input file base name + _summits.fasta',
            'requires_graph': True,
            'arguments':
                [
                    ('peaks_fasta_file', 'Fasta file containing graph peaks.'),
                    ('q_values_base_name', 'Base file name of q values from running '
                                           'the peak caller. Will be peak caller output base name + _qvalues'),
                    ('window_size', 'Optional. Number of basepairs to each side '
                                    'from summit to include. Default is 60.')
                ],
            'method': get_summits
        },
    'project_vg_alignments':
        {
            'help': '',
            'requires_graph': True,
            'arguments':
                [
                    ('alignments_json_file_name', ''),
                    ('linear_path_file_name', ''),
                    ('chromosome', ''),
                    ('out_file_name', 'Writes bed to this file')
                ],
            'method': project_vg_alignments
        },
    'check_pruned_graphs_stats':
        {
            'help': 'Uses vg to check stats of pruned graph vs graphs in graph dir.',
            'arguments':
                [
                    ('graphs_dir', ''),
                    ('chromosomes', 'Comma separated list of chromosomes.'
                                    'Will look for graphs on form graph_dir/[chromosome].pruned.vg')
                ],
            'method': check_pruned_graphs_stats
        },
    'vg_json_alignments_to_intervals':
        {
            'help': 'Reads vg json alignments and converts to interval collection',
            'requires_graph': True,
            'arguments':
                [
                    ('vg_json_file_name', 'Vg json file name'),
                    ('out_file_name', 'Out file name')
                ],
            'method': vg_json_alignments_to_intervals
        },
    'get_intersecting_intervals':
        {
            'help': "Finds intervals in first file that interesects any intervals in second file",
            'requires_graph': True,
            'arguments':
                [
                    ('file1', 'Intervalcollection file 1'),
                    ('file2', 'Intervalcollection file 2'),
                    ('out_file_name', 'File to write intersecting intervals to')
                ],
            'method': get_intersecting_intervals
        }
}


class GraphAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            new_values = obg.Graph.from_numpy_file(values)
        except FileNotFoundError:
            raise GraphNotFoundException()


        setattr(namespace, self.dest, new_values)
        setattr(namespace, "graph_file_name", values)
        try:
            sequencegraph = obg.SequenceGraph.from_file(values + ".sequences")
            setattr(namespace, "sequence_graph", sequencegraph)
            logging.info("Using sequencegraph %s" % (values + ".sequences"))
        except FileNotFoundError:
            logging.info("No sequencegraph found. Will not use sequencegraph.")
            setattr(namespace, "sequence_graph", None)


def run_argument_parser(args):
    # Create parser
    parser = argparse.ArgumentParser(
        description='Graph peak caller',
        prog='graph_peak_caller',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    subparsers = parser.add_subparsers(help='Subcommands')

    for command in interface:

        example = ""
        if "example_run" in interface[command]:
            example = "\nExample: " + interface[command]["example_run"]

        subparser = subparsers.add_parser(
            command,
            help=interface[command]["help"] + example)

        subparser.add_argument('-v', '--verbose',
                               help='Verbosity level. 0 (show only warnings and errors), 1 (show info, '
                                    'warnings and errors), 2 (show everything). Default is 1.', dest='verbose',
                               required=False, metavar='N', const=1, type=int, nargs='?', default=1)

        if 'requires_graph' in interface[command]:
            subparser.add_argument('-g', '--graph', action=GraphAction,
                                   help='Graph file name', dest='graph',
                                   required=True, metavar='file.nobg')

        for argument, help in interface[command]["arguments"]:

            if "/" in argument:
                c = argument.split("/")
                short_command = c[0].strip()
                long_command = c[1].strip()
                assert long_command.startswith("--"), "Long command for %s must start with --" % argument
                assert short_command.startswith("-"), "Short command for %s must start with -" % argument

                required = True
                if "Optional" in help or "optional" in help:
                    required = False

                nargs = None
                if "wildcard" in help:
                    nargs = "+"

                subparser.add_argument(short_command, long_command,
                                       dest=long_command.replace("--", ""), help=help,
                                       required=required, nargs=nargs)
            else:
                subparser.add_argument(argument, help=help)
        subparser.set_defaults(func=interface[command]["method"])


    if len(args) == 0:
        parser.print_help()
        sys.exit(1)


    try:
        args = parser.parse_args(args)
    except GraphNotFoundException:
        logging.error("Specified graph file was not found. Aborting.")
        sys.exit(1)

    set_logging_config(args.verbose)

    if hasattr(args, 'func'):
        try:
            args.func(args)
        except (pyvg.vgobjects.IntervalNotInGraphException, InvalidPileupInterval) as e:
            logging.debug(e)
            logging.error("Found an alignment not compatible with the graph that was used."
                          " Are you sure alignments/intervals are mapped to the same graph that was used?"
                          " Turn on debuging with --verbose 2 to see full log.")
    else:
        parser.help()

if __name__ == "__main__":
    run_argument_parser(sys.argv[1:])
