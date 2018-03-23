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
    run_callpeaks_whole_genome_from_p_values

from graph_peak_caller.analysis.analysis_interface import analyse_peaks_whole_genome,\
    analyse_peaks, peaks_to_fasta, linear_peaks_to_fasta,\
    analyse_manually_classified_peaks, differential_expression,\
    plot_motif_enrichment, get_summits, peaks_to_linear, move_linear_reads_to_graph,\
    find_linear_path, concatenate_sequence_files

from graph_peak_caller.preprocess_interface import \
    count_unique_reads_interface, create_ob_graph,\
    create_linear_map_interface,\
    split_vg_json_reads_into_chromosomes, shift_estimation


def main():
    run_argument_parser(sys.argv[1:])


def version(args):
    print("Graph Peak Caller v1.1.0")


interface = \
{
    'callpeaks':
        {
            'help': 'Call peaks on a single graph.',
            'requires_graph': True,
            'arguments':
                [
                    ('-m/--linear_map', "Optional. Linear map file name. Will look for '"
                                        " files matching graph name if not set or create if"
                                        " no file is found."),
                    ('-s/--sample', 'File name to a vg JSON file or intervalcollection file.'),
                    ('-c/--control', '(Optional) File name to a vg JSON file or intervalcollection file. '
                                     'Only include if a separate control is used.'),
                    ('-n/--out_name', 'Optional. Will be prepended to all output files. Default is nothing.'),
                    ('-f/--fragment_length', 'The fragment length used in this ChIP-seq experiment. If unknown, set to an '
                                        'arbitrary number, e.g. 200. However, for good results, this number should be accurate.'),
                    ('-r/--read_length', 'The read length.'),
                    ('-q/--q_value_threshold', 'Optional. Q value threshold. Default 0.05.')
                ],
            'method': run_callpeaks_interface
        },
    'callpeaks_whole_genome':
        {
            'help': 'Callpeaks on whole genome, using one graph for each chromosome.',
            'arguments':
                [
                    ('chromosomes', 'Comma-separated list of chromosomes to use, e.g. 1,2,X,8,Y'),
                    ('-d/--data_dir', 'Path to data directory containing '
                                      'ob graphs and linear maps.'),
                    ('-s/--sample', 'Sample reads base name. Will use '
                                    'files *_[chromosome].json where * is the base name'),
                    ('-c/--control', 'Optional. Control reads base ame. Will use '
                                     'files *_[chromosome].json where * is the base name'),

                    ('-f/--fragment_length', 'Fragment length.'),
                    ('-r/--read_length', 'Read length'),
                    ('-p/--stop_after_p_values', 'Optional. True or False (default). Whether to '
                                                 'only run until p-value track is '
                                                 'computed (before peak calling)'),
                    ('-u/--unique_reads', 'Number of unique reads. '
                                          'Found by calling count_unique_reads'),
                    ('-g/--genome_size', 'Number of base pairs covered by '
                                         'graphs in total (on a linear genome)'),
                    ('-n/--out_name', 'Optional. Out base name. Prepended to output files.'),
                    ('-D/--keep_duplicates', 'Optional. Set to True in order to keep '
                                             'duplicate input alignments.')
                ],
            'method': run_callpeaks_whole_genome
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
                    ('-r/--read_length', '')
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
                    ('plot_title', 'Title above plot')
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
                    ('graphs_location', 'Will use the graphs *_[chromosome]'),
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
    'peaks_to_linear':
        {
            'help': 'Converts graph peaks to linear peaks.',
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

                subparser.add_argument(short_command, long_command,
                                       dest=long_command.replace("--", ""), help=help,
                                       required=required)
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
