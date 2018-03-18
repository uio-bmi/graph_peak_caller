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

    def pvalues(self, data):
        data.to_sparse_files(
            self._base_name + "pvalues")

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

    def fragment_pileup(self, data):
        data.to_sparse_files(self._base_name+"fragment_pileup")

    def background_track(self, data):
        data.to_sparse_files(self._base_name+"background_track")

    def thersholded(self, data):
        data.to_sparse_files(self._base_name+"thresholded")

    def direct_pileup(self, data):
        data.to_sparse_files(self._base_name+"direct_pileup")

    def touched_nodes(self, data):
        np.save(self._base_name + "touched_nodes.npy",
                np.array(list(data), dtype="int"))

    def add(self, name, data):
        if hasattr(self, name):
            getattr(self, name)(data)
            logging.info("Wrote %s to file", name)
        else:
            logging.info("Skipping reporting of %s", name)

    def get_sub_reporter(self, name):
        return self.__class__(self._base_name + name + "_")
