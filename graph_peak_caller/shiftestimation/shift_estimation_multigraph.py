from ..linear_filter import LinearFilter
from .shiftestimation import Opt, Treatment, PeakModel
import logging


class MultiGraphShiftEstimator(object):
    def __init__(self, position_tuples, genome_size):
        self._position_tuples = position_tuples
        self.genome_size = genome_size

    def get_estimates(self):
        opt = Opt()
        opt.gsize = self.genome_size
        treatment = Treatment(self._position_tuples)
        peakmodel = PeakModel(opt, treatment)
        peakmodel.build()
        return round(peakmodel.d)

    @classmethod
    def from_files(cls, chrom_names, graph_file_names, interval_json_file_names):

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
