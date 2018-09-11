from graph_peak_caller.peakcollection import NonGraphPeakCollection, PeakCollection


def macs_to_graph_peaks(folder):
    for chrom in ["1", "2", "3", "4", "5"]:
        path = NumpyIndexedInterval.from_file("/data/bioinf/tair2/" + chrom  + "_linear_pathv2.interval")
        graph = Graph.from_file("/data/bioinf/tair2/" + chrom + ".nobg")
        macs_peaks = NonGraphPeakCollection.from_fasta(folder + "/macs_sequences_chr%s_summits_unique.fasta" % chrom)
        macs_peaks = PeakCollection.create_from_nongraph_peak_collection(graph, macs_peaks, path)
        macs_peaks.to_file(folder + "/macs_unique_graph_summits_chr%s.intervalcollection" % chrom)

if __name__ == "__main__":
    folder = sys.argv[1]
