import numpy as np
import logging
import offsetbasedgraph as obg
from pyvg.sequences import SequenceRetriever
from . import CallPeaksFromQvalues
from .peakscomparer import PeaksComparerV2, AnalysisResults
from .manually_classified_peaks import CheckOverlapWithManuallyClassifiedPeaks
from .analyse_peaks import LinearRegion
from .analysis.differentialbinding import main
from .analysis.fimowrapper import FimoFile
from .peakcollection import PeakCollection
from .nongraphpeaks import NonGraphPeakCollection


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
    graph = obg.GraphWithReversals.from_numpy_file(args.ob_graph_file_name)

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

    test_name = args.test_name
    fimo_file_name = "fimo_%s_sequences/fimo.txt" % test_name
    peaks_file_name = "%s_max_paths.intervalcollection" % test_name
    subgraphs_file_name = "%s_sub_graphs.graphs.npz" % test_name
    node_ids_file_name = "%s_sub_graphs.nodeids.npz" % test_name
    graph = obg.GraphWithReversals.from_numpy_file(args.graph_name)
    res = main(
        FimoFile.from_file(fimo_file_name),
        PeakCollection.from_file(peaks_file_name, True),
        np.load(subgraphs_file_name),
        np.load(node_ids_file_name),
        graph)
    retriever = SequenceRetriever.from_vg_json_graph(args.vg_graph_name)
    out_f = open(test_name + "_diffexpr.fasta", "w")
    for expr_diff in res:
        main_seq = retriever.get_interval_sequence(expr_diff.main_path)
        out_f.write("> %s %s\n" % (expr_diff.peak_id, expr_diff.main_count))
        out_f.write(main_seq + "\n")
        var_seq = retriever.get_interval_sequence(expr_diff.var_path)
        out_f.write("> %sAlt %s\n" % (expr_diff.peak_id, expr_diff.var_count))
        out_f.write(var_seq + "\n")
    out_f.close()

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


def intervals_to_fasta(args):
    logging.info("Getting sequence retriever")
    retriever = SequenceRetriever.from_vg_graph(args.vg_graph_file_name)
    logging.info("Getting intervals")
    intervals = obg.IntervalCollection.create_generator_from_file(
        args.intervals_file_name)
    logging.info("Writing to fasta")
    CallPeaksFromQvalues.intervals_to_fasta_file(
        intervals, args.out_file_name, retriever)


def linear_peaks_to_fasta(args):
    collection = NonGraphPeakCollection.from_bed_file(
        args.linear_reads_file_name)
    collection.set_peak_sequences_using_fasta(
        fasta_file_location=args.fasta_file)
    collection.save_to_sorted_fasta(args.out_file_name)
    logging.info("Saved sequences to %s" % args.out_file_name)
