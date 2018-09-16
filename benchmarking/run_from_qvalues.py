import sys
import logging
logging.basicConfig(level=logging.DEBUG)
from offsetbasedgraph import Graph, NumpyIndexedInterval
from offsetbasedgraph.vcfmap import load_variant_maps
from graph_peak_caller.postprocess.maxpaths import SparseMaxPaths
from graph_peak_caller.sparsediffs import SparseValues
from graph_peak_caller.peakcollection import PeakCollection

chrom = sys.argv[1]
fragment_length = int(sys.argv[2])

ref = NumpyIndexedInterval.from_file("/data/bioinf/tair2/" + chrom + "_linear_pathv2.interval")


graph = Graph.from_file("/data/bioinf/tair2/" + chrom + ".nobg")
direct = SparseValues.from_sparse_files(chrom + "_direct_pileup")
filtered_peaks = SparseValues.from_sparse_files(chrom + "_hole_cleaned")
variant_map = load_variant_maps(chrom, "/data/bioinf/tair2/")

max_paths, sub_graphs = SparseMaxPaths(filtered_peaks, graph, direct, ref, variant_map).run()
long_maxpaths = [path for path in max_paths if path.length() >= fragment_length]

for max_path in long_max_paths:
    assert max_path.length() > 0, "Max path %s has negative length" % max_path
    score = np.max(self.q_values.get_interval_values(max_path))
    max_path.set_score(score)
    assert not np.isnan(score), "Score %s is nan" % score


PeakCollection(long_maxpaths).to_file(chrom + "_max_paths.intervalcollection", text_file=True)

from graph_peak_caller.peakfasta import PeakFasta
from offsetbasedgraph import SequenceGraph
seqgraph = SequenceGraph.from_file("/data/bioinf/tair2/" + chrom + ".nobg.sequences")
PeakFasta(seqgraph).write_max_path_sequences(chrom + "_sequences.fasta", long_maxpaths)

