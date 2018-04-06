from ..linear_filter import LinearFilter
from .shiftestimation import Opt, Treatment, PeakModel
import logging


class MultiGraphShiftEstimator(object):
    def __init__(self, position_tuples, genome_size,
                 min_fold_enrichment=5, max_fold_enrichment=50):
        self._position_tuples = position_tuples
        self.genome_size = genome_size
        self.min_fold_enrichment = min_fold_enrichment
        self.max_fold_enrichment = max_fold_enrichment

    def get_estimates(self):
        opt = Opt(self.min_fold_enrichment, self.max_fold_enrichment)
        opt.gsize = self.genome_size
        treatment = Treatment(self._position_tuples)
        peakmodel = PeakModel(opt, treatment)
        peakmodel.build()
        return round(peakmodel.d)

    @classmethod
    def from_files(cls, graph_file_names, interval_json_file_names, min_fold_enrichment=5,
                   max_fold_enrichment=50):
        start_positions = {
            "+": {str(i): [] for i, _ in enumerate(graph_file_names)},
            "-": {str(i): [] for i, _ in enumerate(graph_file_names)}
        }
        genome_size = 0
        i = 0
        for graph, intervals in zip(graph_file_names,
                                    interval_json_file_names):
            linear_filter = LinearFilter.from_vg_json_reads_and_graph(
                intervals, graph)

            positions = linear_filter.find_start_positions()
            start_positions["+"][str(i)] = positions["+"]
            start_positions["-"][str(i)] = positions["-"]
            i += 1

            genome_size += linear_filter._indexed_interval.length()

        logging.info("Using genome size %d" % genome_size)
        return cls(start_positions, genome_size, min_fold_enrichment, max_fold_enrichment)

    @classmethod
    def _from_files(cls, graph_file_names, interval_json_file_names):
        start_positions = {
            "+": {name: [] for name in chrom_names},
            "-": {name: [] for name in chrom_names}
        }
        genome_size = 0

        for name, graph, intervals in zip(chrom_names,
                                          graph_file_names,
                                          interval_json_file_names):

            linear_filter = LinearFilter.from_vg_json_reads_and_graph(
                intervals, graph)

            positions = linear_filter.find_start_positions()
            start_positions["+"][name] = positions["+"]
            start_positions["-"][name] = positions["-"]

            genome_size += linear_filter._indexed_interval.length()

        logging.info("Using genome size %d" % genome_size)
        return cls(start_positions, genome_size)
