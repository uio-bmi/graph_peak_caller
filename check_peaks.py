import pysam
from pyvg.alignmentcollection import AlignmentCollection
from offsetbasedgraph import Graph
from graph_peak_caller.peakcollection import PeakCollection
from offsetbasedgraph import NumpyIndexedInterval
import logging

logging.basicConfig(level=logging.DEBUG)

class PeakAnalyser:
    def __init__(self, tf_experiment_dir, data_dir):
        self.experiment_dir = tf_experiment_dir
        self.data_dir = data_dir
        self.bam_file = pysam.AlignmentFile(self.experiment_dir + "/linear_alignments.bam", "rb")
        self.linear_path = NumpyIndexedInterval.from_file(self.data_dir + "/5_linear_pathv2.interval")
        self.graph = Graph.from_file(self.data_dir + "/5.nobg")
        self.alignment_collection = AlignmentCollection.from_file(self.experiment_dir + "/5_alignments.pickle", self.graph)
        self.check_peaks()

    def check_peaks(self):
        peaks = PeakCollection.from_file(self.experiment_dir + "/not_matching_set1.intervals", self.graph) 
        
        i = 0
        for peak in peaks:
            peak.graph = self.graph
            alignments = set([a.strip() for a in self.alignment_collection.get_alignments_on_interval(peak).keys()])
            linear_peak = peak.to_linear_offsets2(self.linear_path)
            linear_alignments = set([a.qname for a in self.bam_file.fetch("5", linear_peak[0], linear_peak[1])])
            i += 1
            print(" ==== Peak %d === " % i)
            print(peak)
            print("Linear: %s" % str(linear_peak))
            print("%d alignments" % len(alignments))
            print("%d linear alignments" % len(linear_alignments))
            print("%d alignments in common" % len(alignments.intersection(linear_alignments)))
            print("%d linear not in graph" % len(linear_alignments.difference(alignments)))
            print("%d graph not in linear" % len(alignments.difference(linear_alignments)))
            print("Not aligned by linear:")
            print(alignments.difference(linear_alignments))


if __name__ == "__main__":
    PeakAnalyser("/var/lib/jenkins/workspace/graph_peak_caller_benchmarks/graph_peak_caller/benchmarking/data/ARABIDOPSIS_ERF115_SRR931836/1/", "/data/bioinf/tair2/")



        




