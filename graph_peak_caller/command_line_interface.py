#!/usr/bin/python3
import argparse
import logging
import sys
import matplotlib as mpl
from graph_peak_caller.densepileup import DensePileup
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

from graph_peak_caller.preprocess_interface import \
    count_unique_reads_interface, create_ob_graph,\
    create_linear_map_interface,\
    split_vg_json_reads_into_chromosomes, shift_estimation


logging.basicConfig(
    stream=sys.stdout, level=logging.WARNING,
    format="%(asctime)s, %(levelname)s: %(message)s")


def main():
    run_argument_parser(sys.argv[1:])


def version(args):
    print("Graph Peak Caller v1.0.4")


def concatenate_sequence_files(args):
    chromosomes = args.chromosomes.split(",")
    out_file_name = args.out_file_name

    if args.is_summits == "True":
        file_endings = "_sequences_summits.fasta"
    else:
        file_endings = "_sequences.fasta"

    all_fasta_entries = []
    for chromosome in chromosomes:
        print("Processing chromosome %s" % chromosome)
        fasta_file = open(chromosome + file_endings)
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
    graph = args.graph
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


def get_summits(args):
    graph = args.graph
    qvalues = DensePileup.from_sparse_files(graph, args.q_values_base_name)
    logging.info("Q values fetched")
    peaks = PeakCollection.from_fasta_file(args.peaks_fasta_file, graph)
    peaks.cut_around_summit(qvalues)
    peaks.to_fasta_file(args.peaks_fasta_file.split(".")[0] + "_summits.fasta", args.sequence_graph)




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
                    ('-n/--out_name', 'Optional. Out base name. Prepended to output files.')
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
    'intervals_to_fasta':
        {
            'help': 'Get sequences for intervals in interval file.',
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
            'help': 'Converts bed file of peaks to fasta',
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
                                        'If False or not set, will search for *_sequences.fasta')
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
                    ('chromosomes', 'Graph base names. Set to empty string if only single  graph is being used. If whole-genome, use comma-separated list of chromosomes to use, e.g. 1,2,X,8,Y'),
                    ('ob_graphs_location', 'Location of graph files'),
                    ('sample_reads_base_name', 'Will use files [base_name][chromosome].json'),
                ],
            'method': shift_estimation
        },
    'get_summits':
        {
            'help': 'Get summit around peaks.',
            'requires_graph': True,
            'arguments':
                [
                    ('peaks_fasta_file', 'Fasta file containing graph peaks.'),
                    ('q_values_base_name', 'Base file name of q values')
                ],
            'method': get_summits
        }
}



class GraphAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        new_values = obg.Graph.from_numpy_file(values)

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

        if 'requires_graph' in interface[command]:
            subparser.add_argument('-g', '--graph', action=GraphAction,
                                   help='Graph file name', dest='graph',
                                   required=True, metavar='GRAPH')

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

    args = parser.parse_args(args)
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.help()

if __name__ == "__main__":
    run_argument_parser(sys.argv[1:])
