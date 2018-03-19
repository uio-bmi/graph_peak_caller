import logging
from .sparsegraphpileup import SamplePileupGenerator


def get_fragment_pileup(graph, input_intervals, info, reporter=None):
    logging.info("Creating fragment pileup")
    spg = SamplePileupGenerator(graph, info.fragment_length-info.read_length)
    return spg.run(input_intervals, reporter)
