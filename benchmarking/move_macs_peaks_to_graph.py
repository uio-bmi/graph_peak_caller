from graph_peak_caller.analysis.nongraphpeaks import NonGraphPeakCollection
from graph_peak_caller.peakcollection import PeakCollection

import sys
from offsetbasedgraph import Graph, NumpyIndexedInterval


for chrom in sys.argv[1].split(","):
    linear_path = NumpyIndexedInterval.from_file("/data/bioinf/tair2/" + chrom + "_linear_pathv2.interval")
    graph = Graph.from_file("/data/bioinf/tair2/" + chrom + ".nobg")
    print("Chrom " + chrom)
    peaks = NonGraphPeakCollection.from_bed_file("macs_peaks_chr" + chrom + ".bed", 60) 
    print(len(peaks.peaks))
    graph_peaks = PeakCollection.create_from_nongraph_peak_collection(graph, peaks, linear_path)
    graph_peaks.to_file(chrom + "_macs_all_summits.intervalcollection", text_file=True)
    print("Wrote to " + chrom + "_macs_all_summits.intervalcollection")


    






