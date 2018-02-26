
from .linear_filter import LinearFilter
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

    def to_linear_bed_file(self, file_name, read_length):
        f = open(file_name, "w")
        for direction in self._position_tuples.keys():
            for chrom, positions in self._position_tuples[direction].items():
                for pos in positions:
                    #print("Writing position %d" % pos)

                    start = pos
                    end = pos + read_length

                    if direction == "-":
                        start = pos - read_length
                        end = pos

                    f.writelines(["chr%s\t%d\t%d\t.\t0\t%s\n" % \
                                  (chrom, start, end, direction)])
        f.close()


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
