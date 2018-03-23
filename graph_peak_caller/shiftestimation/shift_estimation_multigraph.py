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
