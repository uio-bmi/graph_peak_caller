import logging
import numpy as np

from .peakcollection import PeakCollection


class Reporter:
    def __init__(self, base_name):
        self._base_name = base_name

    def sub_graphs(self, data):
        np.savez(self._base_name + "sub_graphs.graphs",
                 **{"peak%s" % i: p._graph for i, p in enumerate(data)})
        np.savez(self._base_name + "sub_graphs.nodeids",
                 **{"peak%s" % i: p._node_ids for i, p in enumerate(data)})

    def qvalues(self, data):
        data.to_sparse_files(
            self._base_name + "qvalues")

    def all_max_paths(self, data):
        PeakCollection(data).to_file(
            self._base_name+"all_max_paths.intervalcollection",
            text_file=True)

    def max_paths(self, data):
        PeakCollection(data).to_file(
            self._base_name+"max_paths.intervalcollection",
            text_file=True)

    def hole_cleaned(self, data):
        data.to_sparse_files(self._base_name+"hole_cleaned")

    def thersholded(self, data):
        data.to_sparse_files(self._base_name+"thresholded")

    def add(self, name, data):
        if hasattr(self, name):
            getattr(self, name)(data)
            logging.info("Wrote %s to file", name)
        else:
            logging.info("Skipping reporting of %s", name)
