import logging
from .sparsegraphpileup import SamplePileupGenerator


def get_fragment_pileup(graph, input_intervals, info, reporter=None, ):
    logging.info("Creating fragment pileup, using fragment length %d "
                 "and read length %d" % (info.fragment_length, info.read_length))
    spg = SamplePileupGenerator(graph, info.fragment_length-info.read_length)
    return spg.run(input_intervals, reporter)

def get_fragment_pileup_atac(graph, input_intervals, info, reporter=None):
    logging.info("Creating fragment pileup, using fragment length %d "
                 "and read length %d" % (info.fragment_length, info.read_length))
    spg = SamplePileupGenerator(graph, info.fragment_length//2-info.read_length+info.fragment_length % 2, info.fragment_length//2)
    return spg.run(input_intervals, reporter)

