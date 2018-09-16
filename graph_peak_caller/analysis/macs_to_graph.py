from graph_peak_caller.peakcollection import NonGraphPeakCollection, PeakCollection
from offsetbasedgraph import NumpyIndexedInterval, Graph

def macs_to_graph_peaks(folder):
    for chrom in ["1", "2", "3", "4", "5"]:
        path = NumpyIndexedInterval.from_file("/data/bioinf/tair2/" + chrom  + "_linear_pathv2.interval")
        graph = Graph.from_file("/data/bioinf/tair2/" + chrom + ".nobg")
        macs_peaks = PeakCollection.from_fasta_file(folder + "/macs_sequences_chr%s_summits_unique.fasta" % chrom, graph)
        macs_peaks.to_file(folder + "/%s_macs_unique_graph_summits.intervalcollection" % chrom, True)

if __name__ == "__main__":
    import sys
    folder = sys.argv[1]
    macs_to_graph_peaks(folder)
