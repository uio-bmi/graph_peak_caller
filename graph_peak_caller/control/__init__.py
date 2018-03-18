import logging

from .controlgenerator import SparseControl
from .linearmap import LinearMap


logger = logging.getLogger("background")


def get_background_track(graph, linear_map_name, extensions, fragment_length,
                         intervals, touched_nodes=None, global_min=None):
    sc = SparseControl(linear_map_name, graph, extensions,
                       fragment_length, touched_nodes)
    if global_min is not None:
        sc.set_min_value(global_min)
    return sc.create(intervals)


def get_background_track_from_control(
        graph, linear_map_name, fragment_length,
        intervals, touched_nodes=None, global_min=None):
    logger.info("Creating background from control reads")
    extensions = [fragment_length, 1000, 10000]
    sc = SparseControl(linear_map_name, graph, extensions,
                       fragment_length, touched_nodes)
    if global_min is not None:
        sc.set_min_value(global_min)
    return sc.create(intervals)


def get_background_track_from_input(
        graph, linear_map_name, fragment_length,
        intervals, touched_nodes=None, global_min=None):
    logger.info("Creating background from input reads")
    extensions = [10000]
    sc = SparseControl(linear_map_name, graph, extensions,
                       fragment_length, touched_nodes)
    if global_min is not None:
        sc.set_min_value(global_min)
    return sc.create(intervals)


def scale_tracks(fragment_pileup, background_track, ratio):
    logger.info("Scaling tracks to ratio: %d" % ratio)
    if ratio == 1:
        return

    if ratio > 1:
        logging.warning("More reads in sample than in control")
        fragment_pileup *= 1/ratio
    else:
        logging.info("Scaling control pileup down with ratio %.3f" % ratio)
        background_track *= ratio
