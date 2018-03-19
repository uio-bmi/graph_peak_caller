import logging

from .controlgenerator import SparseControl
from .linearmap import LinearMap


def get_background_track(graph, intervals, config, extensions,
                         touched_nodes=None):
    sc = SparseControl(config.linear_map_name, graph, extensions,
                       config.fragment_length, touched_nodes)
    if config.global_min is not None:
        sc.set_min_value(config.global_min)
    return sc.create(intervals)


def get_background_track_from_control(
        graph, intervals, config, touched_nodes=None):
    logging.info("Creating background from control reads")
    extensions = [config.fragment_length, 1000, 10000]
    return get_background_track(graph, intervals, config,
                                extensions, touched_nodes)


def get_background_track_from_input(
        graph, intervals, config, touched_nodes=None):
    logging.info("Creating background from input reads")
    extensions = [10000]
    return get_background_track(graph, intervals, config,
                                extensions, touched_nodes)


def scale_tracks(fragment_pileup, background_track, ratio):
    logging.info("Scaling tracks to ratio: %d" % ratio)
    if ratio == 1:
        return

    if ratio > 1:
        logging.warning("More reads in sample than in control")
        fragment_pileup *= 1/ratio
    else:
        logging.info("Scaling control pileup down with ratio %.3f" % ratio)
        background_track *= ratio
