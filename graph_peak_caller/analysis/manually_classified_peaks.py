from graph_peak_caller.peakcollection import PeakCollection
import pyvg
import logging
from graph_peak_caller.util import create_linear_path
import offsetbasedgraph as obg

class CheckOverlapWithManuallyClassifiedPeaks:

    def __init_(self, chromosome, reads, true_peaks_file):
        # Filter out reads that are not in the regions in region_bed_file_name
        self.chromosome = chromosome
        self.reads = reads

    @classmethod
    def from_graph_peaks_in_fasta(cls, graph, vg_graph_json_file_name, chromosome, fasta_file_name, regions_bed_file, true_peaks_file):
        reads = PeakCollection.from_fasta_file(fasta_file_name, graph=graph)
        vg_graph = pyvg.vg.Graph.create_from_file(vg_graph_json_file_name, limit_to_chromosomes=chromosome)
        logging.info("Finding linear path")

        linear_path_file = "linear_path_%s.intervalcollection" % chromosome
        try:
            linear_path = obg.IntervalCollection.from_file(linear_path_file, text_file=True).intervals[0]
            linear_path = linear_path.to_indexed_interval()
        except FileNotFoundError:
            linear_path = create_linear_path(graph, vg_graph, path_name=chromosome, write_to_file=linear_path_file)

        linear_path.graph = graph

        filtered_reads = []

        # Convert regions to intervals in graph
        logging.info("Converting regions to regions in graph")
        graph_regions = []
        bed_file = open(regions_bed_file)
        for line in bed_file:
            print(line)
            line = line.split()
            chr = line[0]
            start = int(line[1])
            end = int(line[2])

            if chr != "chr%s" % chromosome:
                logging.info("Skipping %s, %d, %d" % (chr, start, end))
                continue

            graph_interval = linear_path.get_subinterval(start, end)
            graph_regions.append(graph_interval)

        assert len(graph_regions) > 0, " Found not graph regions for chr %d" % chromosome
        graph_regions = PeakCollection(graph_regions)

        # Filter out reads not overlapping with linear regions
        for read in reads:
            n_overlapping = graph_regions.get_overlapping_intervals(read, minimum_overlap=1)

            if n_overlapping:
                filtered_reads.append(read)

        logging.info("Found %d reads in graph regions" % len(filtered_reads))

        return cls(chromosome, reads, true_peaks_file)
