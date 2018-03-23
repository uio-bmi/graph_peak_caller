import numpy as np
import logging
import offsetbasedgraph as obg
import pyvg
from pyvg.sequences import SequenceRetriever
from .peakscomparer import PeaksComparerV2, AnalysisResults
from .manually_classified_peaks import \
    CheckOverlapWithManuallyClassifiedPeaks
from .analyse_peaks import LinearRegion
from .differentialbinding import main
from .fimowrapper import FimoFile
from ..peakcollection import Peak, PeakCollection
from .nongraphpeaks import NonGraphPeakCollection
from .motifenrichment import plot_true_positives
from ..peakfasta import PeakFasta
from ..mindense import DensePileup
from ..sparsediffs import SparseValues
from .util import create_linear_path


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

        if args.use_input_file_pattern is not None:
            assert "[chrom]" in args.use_input_file_pattern, \
                "Use [chrom] to specify where chromosome should be replaced."
            try:
                fasta_file = open(args.use_input_file_pattern.replace("[chrom]", chromosome))
            except FileNotFoundError:
                logging.info("Did not find file matching patterh %s. Does those files exist?" % args.use_input_file_pattern)
                raise
        else:
            logging.info("Guessing file name since use_input_file_pattern is not specified")
            fasta_file = open(chromosome + file_endings)

        for line in fasta_file:
            if line.startswith(">"):
                all_fasta_entries.append([line, None])
            else:
                # This is sequence, add to prev entry
                all_fasta_entries[-1][1] = line
        fasta_file.close()

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
    qvalues = SparseValues.from_sparse_files(args.q_values_base_name)
    dense_qvalues = qvalues.to_dense_pileup(size=graph.number_of_basepairs())
    qvalues = DensePileup(graph, dense_qvalues)
    logging.info("Q values fetched")
    peaks = PeakCollection.from_fasta_file(args.peaks_fasta_file, graph)

    if args.window_size is not None:
        window = int(args.window_size)
    else:
        window = 60

    peaks.cut_around_summit(qvalues, n_base_pairs_around=window)
    out_file_name = args.peaks_fasta_file.split(".")[0] + "_summits.fasta"
    peaks.to_fasta_file(out_file_name, args.sequence_graph)
    logging.info("Wrote summits to " + out_file_name)


def analyse_peaks_whole_genome(args):
    chromosomes = args.chromosomes.split(",")
    results = AnalysisResults()
    for chrom in chromosomes:
        graph = obg.GraphWithReversals.from_numpy_file(
            args.graphs_dir + chrom + ".nobg")
        logging.info("Reading linear path")
        linear_path = obg.NumpyIndexedInterval.from_file(
            args.graphs_dir + chrom + "_linear_pathv2.interval")

        analyser = PeaksComparerV2(
            graph,
            args.results_dir + "macs_sequences_chr%s_summits.fasta" % chrom,
            args.results_dir + "%s_sequences_summits.fasta" % chrom,
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
    logging.info("Wrote results as pickle to %s and as text to %s"
                 % (args.out_file + ".pickle", args.out_file + ".txt"))


def analyse_manually_classified_peaks(args):
    for chromosome in args.chromosomes.split():
        graph_file_name = args.graphs_location + "/" + chromosome + ".nobg"
        graph = obg.GraphWithReversals.from_numpy_file(graph_file_name)
        CheckOverlapWithManuallyClassifiedPeaks.from_graph_peaks_in_fasta(
                graph,
                args.graphs_location + "/" + chromosome + ".json",
                chromosome,
                args.reads_base_name + chromosome + "_sequences.fasta",
                args.regions_file,
                args.manually_classified_peaks_file)


def analyse_peaks(args):
    graph = args.graph

    end = int(args.graph_end)
    if end == 0:
        region = None
    else:
        region = LinearRegion(args.graph_chromosome,
                              int(args.graph_start), end)
    linear_path = obg.NumpyIndexedInterval.from_file(
        args.linear_path_name)
    PeaksComparerV2(
        graph,
        args.linear_peaks_fasta_file_name,
        args.graph_peaks_fasta_file_name,
        args.linear_peaks_fimo_results_file,
        args.graph_peaks_fimo_results_file,
        linear_path,
        region=region)


def differential_expression(args):
    logging.info("Running differential expression.")
    test_name = args.test_name
    fimo_file_name = args.fimo_file_name
    peaks_file_name = "%s_max_paths.intervalcollection" % test_name
    subgraphs_file_name = "%s_sub_graphs.graphs.npz" % test_name
    node_ids_file_name = "%s_sub_graphs.nodeids.npz" % test_name
    graph = args.graph

    subgraphs = np.load(subgraphs_file_name)
    logging.info("Found %d subgraphs" % len(subgraphs.keys()))
    node_ids = np.load(node_ids_file_name)

    res = main(
        FimoFile.from_file(fimo_file_name),
        PeakCollection.from_file(peaks_file_name, True),
        subgraphs,
        node_ids,
        graph)
    retriever = obg.SequenceGraph.from_file(
        args.graph_file_name + ".sequences")
    out_file_name = test_name + "_diffexpr.fasta"
    out_f = open(out_file_name, "w")
    n = 0
    for expr_diff in res:
        n += 1
        main_seq = retriever.get_interval_sequence(expr_diff.main_path)
        out_f.write("> %s %s\n" % (expr_diff.peak_id, expr_diff.main_count))
        out_f.write(main_seq + "\n")
        var_seq = retriever.get_interval_sequence(expr_diff.var_path)
        out_f.write("> %sAlt %s\n" % (expr_diff.peak_id, expr_diff.var_count))
        out_f.write(var_seq + "\n")
    logging.info("Wrote %d lines to %s" % (n, out_file_name))
    out_f.close()


def plot_motif_enrichment(args):
    fasta1 = args.fasta1
    fasta2 = args.fasta2
    meme = args.meme_motif_file

    plot_true_positives(
        [
            ("Graph Peak Caller", fasta1),
            ("MACS2", fasta2)
        ],
        meme,
        plot_title=args.plot_title.replace("ARABIDOPSIS_", ""),
        save_to_file=args.out_figure_file_name,
        run_fimo=args.run_fimo == "True"
    )


def peaks_to_fasta(args):
    logging.info("Getting sequence retriever")
    retriever = obg.SequenceGraph.from_file(args.sequence_graph)
    logging.info("Getting intervals")
    intervals = PeakCollection.create_generator_from_file(
        args.intervals_file_name)
    logging.info("Writing to fasta")
    PeakFasta(retriever).save_intervals(args.out_file_name,
                                        intervals)


def linear_peaks_to_fasta(args):
    collection = NonGraphPeakCollection.from_bed_file(
        args.linear_reads_file_name)
    collection.set_peak_sequences_using_fasta(
        fasta_file_location=args.fasta_file)
    collection.save_to_sorted_fasta(args.out_file_name)
    logging.info("Saved sequences to %s" % args.out_file_name)

    window = 60
    if hasattr(args, "window"):
        if args.window is not None:
            window = int(args.window)
            logging.info("Using window size of %d" % window)

    summits = NonGraphPeakCollection.from_bed_file(
        args.linear_reads_file_name, cut_around_summit=window)
    summits.set_peak_sequences_using_fasta(
        fasta_file_location=args.fasta_file)
    out_name = args.out_file_name.split(".")[0] + "_summits.fasta"
    summits.save_to_sorted_fasta(out_name)
    logging.info("Saved summits to %s" % out_name)
