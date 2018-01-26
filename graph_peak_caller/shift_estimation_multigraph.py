
from pyvg.linear_filter import LinearFilter
from .shiftestimation import Opt, Treatment


def MultiGraphShiftEstimator(object):

    def __init__(self, position_tuples):
        self._position_tuples = position_tuples

    def get_fragment_length(self):
        opt = Opt()
        treatment = Treatment(self._position_tuples)

    @classmethod
    def from_files(chrom_names, graph_file_names, interval_json_file_names):

        start_positions = {
            "+": {name: [] for name in chrom_names},
            "-": {name: [] for name in chrom_names}
        }

        for name, graph, intervals in zip(chrom_names,
                                          graph_file_names,
                                          interval_json_file_names):

            linear_filter = LinearFilter.from_vg_json_reads_and_graph(
                intervals, graph)

            positions = linear_filter.find_start_positions()
            start_positions["+"][name] = positions["+"]
            start_positions["-"][name] = positions["-"]

