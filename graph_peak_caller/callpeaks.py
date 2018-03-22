import logging
import numpy as np
from .mindense import DensePileup
from .sample import get_fragment_pileup
from .control import get_background_track_from_control,\
    get_background_track_from_input, scale_tracks
from .sparsepvalues import PValuesFinder, PToQValuesMapper, QValuesFinder
from .postprocess import HolesCleaner, SparseMaxPaths
from .sparsediffs import SparseValues
import json
import sys

class Configuration:
    def __init__(self):
        self.read_length = None
        self.fragment_length = None
        self.linear_map_name = None
        self.has_control = False
        self.q_values_threshold = 0.05
        self.global_min = None
        self.keep_duplicates = False

    def copy(self):
        o = Configuration()
        o.read_length = self.read_length
        o.fragment_length = self.fragment_length
        o.linar_map_name = self.linear_map_name
        o.has_control = self.has_control
        o.q_values_threshold = self.q_values_threshold
        o.global_min = self.global_min
        o.keep_duplicates = self.keep_duplicates
        return o


class CallPeaks(object):

    def __init__(self, graph, config, reporter):
        self.graph = graph
        self.config = config
        assert self.config.fragment_length > self.config.read_length, \
            "Read length %d is larger than fragment size" % (self.config.read_length)
        self.info = config
        self._reporter = reporter

    def run_pre_callpeaks(self, input_reads, control_reads):
        try:
            sample_pileup = get_fragment_pileup(
                self.graph, input_reads, self.config,
                self._reporter)
            background_func = get_background_track_from_control
        except json.decoder.JSONDecodeError as e:
            logging.debug(e)
            logging.critical("Error when parsing json reads. Input file is not valid JSON."
                             "Turn on debugging (--verbose 2) for more debug info.")
            sys.exit(1)

        if not self.config.has_control:
            background_func = get_background_track_from_input
        control_pileup = background_func(self.graph, control_reads,
                                         self.config,
                                         sample_pileup.touched_nodes)
        scale_tracks(sample_pileup, control_pileup,
                     input_reads.n_reads/control_reads.n_reads)

        self._reporter.add("fragment_pileup", sample_pileup)
        self._reporter.add("touched_nodes", sample_pileup.touched_nodes)
        self._reporter.add("background_track", control_pileup)
        self.touched_nodes = sample_pileup.touched_nodes
        self.control_pileup = control_pileup
        self.sample_pileup = sample_pileup

    def get_p_values(self):
        assert self.sample_pileup is not None
        assert self.control_pileup is not None
        self.p_values_pileup = PValuesFinder(
            self.sample_pileup, self.control_pileup).get_p_values_pileup()
        self.p_values_pileup.track_size = self.graph.node_indexes[-1]
        self._reporter.add("pvalues", self.p_values_pileup)
        self.sample_pileup = None
        self.control_pileup = None

    def get_p_to_q_values_mapping(self):
        assert self.p_values_pileup is not None
        finder = PToQValuesMapper.from_p_values_pileup(
            self.p_values_pileup)
        self.p_to_q_values_mapping = finder.get_p_to_q_values()

    def get_q_values(self):
        assert self.p_values_pileup is not None
        assert self.p_to_q_values_mapping is not None
        finder = QValuesFinder(
            self.p_values_pileup,
            self.p_to_q_values_mapping)

        self.q_values_pileup = finder.get_q_values()
        self.q_values_pileup.track_size = self.p_values_pileup.track_size
        self._reporter.add("qvalues", self.q_values_pileup)

    def call_peaks_from_q_values(self):
        assert self.q_values_pileup is not None
        caller = CallPeaksFromQvalues(
            self.graph, self.q_values_pileup,
            self.config, self._reporter,
            touched_nodes=self.touched_nodes,
            config=self.config
            )
        caller.callpeaks()
        self.max_path_peaks = caller.max_paths

    def run(self, input_intervals, control_intervals):
        self.run_to_p_values(input_intervals, control_intervals)
        self.get_p_to_q_values_mapping()
        self.get_q_values()
        self.call_peaks_from_q_values()

    def run_to_p_values(self, input_intervals, control_intervals):
        self.run_pre_callpeaks(input_intervals, control_intervals)
        self.get_p_values()


class CallPeaksFromQvalues:
    def __init__(self, graph, q_values_pileup,
                 experiment_info, reporter,
                 cutoff=0.1, raw_pileup=None, touched_nodes=None,
                 config=None, q_values_max_path=False):

        self.graph = graph
        self.q_values = q_values_pileup
        self.info = experiment_info
        self.cutoff = cutoff
        self.raw_pileup = raw_pileup
        self.touched_nodes = touched_nodes
        self.save_tmp_results_to_file = False
        self.graph_is_partially_ordered = False
        self.q_values_max_path = q_values_max_path

        if config is not None:
            self.cutoff = config.q_values_threshold
        self._reporter = reporter
        # self.info.to_file(self.out_file_base_name + "experiment_info.pickle")
        logging.info("Using p value cutoff %.4f" % self.cutoff)

    def __threshold(self):
        threshold = -np.log10(self.cutoff)
        logging.info("Thresholding peaks on q value %.4f" % threshold)
        self.pre_processed_peaks = self.q_values.threshold_copy(threshold)
        self._reporter.add("thresholded", self.pre_processed_peaks)

    def __postprocess(self):
        logging.info("Filling small Holes")
        self.pre_processed_peaks = HolesCleaner(
            self.graph,
            self.pre_processed_peaks,
            self.info.read_length,
            self.touched_nodes
        ).run()
        self._reporter.add("hole_cleaned", self.pre_processed_peaks)

        logging.info("Not removing small peaks")
        self.filtered_peaks = self.pre_processed_peaks

    def __get_max_paths(self):
        logging.info("Getting maxpaths")
        if not self.q_values_max_path:
            file_name = self._reporter._base_name+"direct_pileup"
            logging.info("Reading raw direct pileup from file %s" % file_name)
            _pileup = SparseValues.from_sparse_files(file_name)
            self.raw_pileup = _pileup
        else:
            # raise NotImplementedError()
            _pileup = self.q_values

        assert(self.graph.uses_numpy_backend)
        logging.info("Running Sparse Max Paths")
        max_paths, sub_graphs = SparseMaxPaths(
            self.filtered_peaks, self.graph, _pileup).run()

        self._reporter.add("all_max_paths", max_paths)
        logging.info("All max paths found")

        self.q_values = DensePileup(
            self.graph, self.q_values.to_dense_pileup(
                self.graph.node_indexes[-1]))

        for max_path in max_paths:
            assert max_path.length() > 0, "Max path %s has negative length" % max_path
            score = np.max(self.q_values.get_interval_values(max_path))
            max_path.set_score(score)
            assert not np.isnan(score), "Score %s is nan" % score

        pairs = list(zip(max_paths, sub_graphs))
        pairs.sort(key=lambda p: p[0].score, reverse=True)
        logging.info("N unfiltered peaks: %s", len(max_paths))
        pairs = [p for p in pairs if
                 p[0].length() >= self.info.fragment_length]
        logging.info("N filtered peaks: %s", len(pairs))
        self._reporter.add("sub_graphs", [pair[1] for pair in pairs])
        self.max_paths = [p[0] for p in pairs]
        self._reporter.add("max_paths", self.max_paths)

    def callpeaks(self):
        logging.info("Calling peaks")
        self.__threshold()
        self.__postprocess()
        self.__get_max_paths()

