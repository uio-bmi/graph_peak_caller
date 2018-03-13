import logging
from pyvg.conversion import vg_json_file_to_interval_collection,\
    json_file_to_obg_numpy_graph
import offsetbasedgraph as obg
from graph_peak_caller.shiftestimation.shift_estimation_multigraph import MultiGraphShiftEstimator
from graph_peak_caller.util import create_linear_map
from graph_peak_caller.multiplegraphscallpeaks import MultipleGraphsCallpeaks


def count_unique_reads_interface(args):
    chromosomes = args.chromosomes.split(",")
    graph_file_names = [args.graphs_location + chrom for chrom in chromosomes]
    reads_file_names = [args.reads_base_name + chrom + ".json"
                        for chrom in chromosomes]

    count_unique_reads(chromosomes, graph_file_names,
                       reads_file_names)

def shift_estimation(args):
    chromosomes = args.chromosomes.split(",")
    graphs = [args.ob_graphs_location + "/" + chrom + ".nobg" for chrom in chromosomes]
    print(graphs)
    logging.info("Will try to use graphs %s" % graphs)
    sample_file_names = [args.sample_reads_base_name + chrom + ".json"
                         for chrom in chromosomes]
    logging.info("Will use reads from %s" % sample_file_names)

    estimator = MultiGraphShiftEstimator.from_files(
        chromosomes, graphs, sample_file_names)

    estimator.to_linear_bed_file("linear_bed.bed", read_length=36)

    d = estimator.get_estimates()
    logging.info("Found shift: %d" % d)

def count_unique_reads(chromosomes, graph_file_names, reads_file_names):
    graphs = (obg.GraphWithReversals.from_numpy_file(f)
              for f in graph_file_names)
    reads = (vg_json_file_to_interval_collection(f, graph)
             for f, graph in zip(reads_file_names, graphs))

    unique_reads = MultipleGraphsCallpeaks.count_number_of_unique_reads(reads)
    print(unique_reads)


def create_ob_graph(args):
    logging.info("Creating obgraph")
    ob_graph = json_file_to_obg_numpy_graph(args.vg_json_file_name, 0)
    logging.info("Writing ob graph to file")
    out_name = args.out_file_name if args.out_file_name is not None \
                else args.vg_json_file_name.split(".")[0] + ".nobg"
    ob_graph.to_numpy_file(out_name)

    logging.info("Creating sequence graph")
    from offsetbasedgraph import SequenceGraph
    sequence_graph = SequenceGraph.create_empty_from_ob_graph(ob_graph)
    sequence_graph.set_sequences_using_vg_json_graph(args.vg_json_file_name)
    out_file_name = out_name + ".sequences"
    sequence_graph.to_file(out_file_name)
    logging.info("Wrote sequence graph to %s" % out_file_name)


def create_linear_map_interface(args):
    graph = args.graph
    # graph.convert_to_dict_backend()

    out_name = args.out_file_base_name if args.out_file_base_name is not None \
        else args.graph_file_name.split(".")[0] + "_linear_map.npz"
    create_linear_map(graph, out_name)
    logging.info("Wrote linear map to file %s" % out_name)


def split_vg_json_reads_into_chromosomes(args):
    reads_base_name = args.vg_json_reads_file_name.split(".")[0]
    logging.info("Will write reads to files %s_[chromosome].json",
                 reads_base_name)

    chromosomes = args.chromosomes.split(",")
    chromosome_limits = {}
    logging.info("Found the following chromosome ranges:")
    for chrom in chromosomes:
        start_end = open(
            args.range_files_base_name + "node_range_" + chrom + ".txt")
        start_end = start_end.read().split(":")
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
