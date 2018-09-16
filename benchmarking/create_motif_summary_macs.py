import logging
logging.basicConfig(level=logging.INFO)
from offsetbasedgraph import Graph, NumpyIndexedInterval
from graph_peak_caller.peakcollection import NonGraphPeakCollection, PeakCollection

# Finds the graph peak IDs where the corresponding macs peaks has a motif match

def get_id_mapping(graph_peaks, macs_peaks):
    # Finds mapping from macs peaks to graph peaks where macs peaks overlap graph peak
    logging.info("Getting matching peaks")
    mapping = {}
    visited = set()

    for peak in graph_peaks.intervals:
        touching = macs_peaks.approx_contains_part_of_interval(
            peak, visited)
        if touching:
            visited.add(touching[0].unique_id)
            mapping[touching[0].unique_id] = peak.unique_id

    return mapping

out_file = open("motif_summary_graph_matching_macs.tsv", "w") 
for chrom in ["1", "2", "3", "4", "5"]: 
    logging.info("Chromosome %s" % chrom)
    path = NumpyIndexedInterval.from_file("/data/bioinf/tair2/" + chrom  + "_linear_pathv2.interval")
    graph = Graph.from_file("/data/bioinf/tair2/" + chrom + ".nobg")
    macs_peaks = NonGraphPeakCollection.from_fasta("macs_sequences_chr" + chrom + "_summits.fasta")
    macs_peaks = PeakCollection.create_from_nongraph_peak_collection(graph, macs_peaks, path)
    macs_peaks.create_node_index()
    graph_peaks = PeakCollection.from_fasta_file(chrom + "_sequences_summits.fasta")
    graph_peaks.create_node_index()
    macs_motif_matches = set([line.split("\t")[2] for line in open("fimo_macs_chr" + chrom + "/fimo.txt") if not line.startswith("#")])
    graph_motif_matches = set([line.split("\t")[2] for line in open("fimo_graph_chr" + chrom + "/fimo.txt") if not line.startswith("#")])

    mapping = get_id_mapping(graph_peaks, macs_peaks)

    motif_ids = [mapping[peak.unique_id] for peak in macs_peaks.intervals if 
                peak.unique_id in macs_motif_matches and peak.unique_id in mapping]
    out_file.write("%s\t%s\n" % (chrom, ",".join(motif_ids)))

out_file.close()
logging.info("Wrote motif summary to motif_summary_graph_matching_macs.tsv")
