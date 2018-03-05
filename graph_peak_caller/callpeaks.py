import logging
import numpy as np
from offsetbasedgraph import IntervalCollection, DirectedInterval
from .densepileup import DensePileup
from .areas import BinaryContinousAreas, BCACollection
from .peakscores import ScoredPeak
from .peakcollection import PeakCollection
IntervalCollection.interval_class = DirectedInterval
from .experiment_info import ExperimentInfo
from .subgraphcollection import SubgraphCollectionPartiallyOrderedGraph
from .peakcollection import Peak
# from memory_profiler import profile
from .sampleandcontrolcreator import SampleAndControlCreator
from .directsamplepileup import DirectPileup
from .sparsepvalues import PValuesFinder, PToQValuesMapper, QValuesFinder
from .sparseholecleaner import HolesCleaner
from .sparsediffs import SparseValues
from .sparsemaxpaths import SparseMaxPaths


class Configuration:
    def __init__(self, skip_filter_duplicates=False,
                 graph_is_partially_ordered=False,
                 skip_read_validation=False,
                 save_tmp_results_to_file=False,
                 p_val_cutoff=0.1,
                 use_global_min_value=None):
        self.skip_filter_duplicates = skip_filter_duplicates
        self.graph_is_partially_ordered = graph_is_partially_ordered
        self.skip_read_validation = skip_read_validation
        self.save_tmp_results_to_file = save_tmp_results_to_file
        self.p_val_cutoff = p_val_cutoff
        self.use_global_min_value = use_global_min_value

    @classmethod
    def default(cls):
        return cls()


class CallPeaks(object):

    def __init__(self, graph, out_file_base_name):
        self.graph = graph
        self.out_file_base_name = out_file_base_name

        self.sample_intervals = None
        self.control_intervals = None
        self.sample_pileup = None
        self.control_pileup = None
        self.p_values_pileup = None
        self.p_to_q_values_mapping = None
        self.q_values_pileup = None
        self.peaks_as_subgraphs = None
        self.touched_nodes = None
        self.max_path_peaks = None

    def run_pre_callpeaks(self, has_control, experiment_info,
                          linear_map, configuration=None):
        if configuration is None:
            configuration = Configuration.default()
            logging.warning("Config is not set. Setting to default")
        creator = SampleAndControlCreator(
            self.graph,
            self.sample_intervals,
            self.control_intervals,
            experiment_info,
            out_file_base_name=self.out_file_base_name,
            has_control=has_control,
            linear_map=linear_map,
            configuration=configuration
            )
        creator.run()
        self.sample_pileup = creator._sample_pileup
        self.control_pileup = creator._control_pileup
        self.touched_nodes = creator.touched_nodes

    def get_p_values(self):
        assert self.sample_pileup is not None
        assert self.control_pileup is not None
        self.p_values_pileup = PValuesFinder(
            self.sample_pileup, self.control_pileup).get_p_values_pileup()
        self.p_values_pileup.track_size = self.graph.node_indexes[-1]
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
        finder = QValuesFinder(self.p_values_pileup,
                               self.p_to_q_values_mapping)
        self.q_values_pileup = finder.get_q_values()
        self.q_values_pileup.track_size = self.p_values_pileup.track_size
        self.q_values_pileup.to_sparse_files(
            self.out_file_base_name + "qvalues")
        logging.info("Wrote q values to sparse files with base name %s " % \
                     (self.out_file_base_name + "qvalues"))
        # q_values = DensePileup(self.graph)
        #q_values.data._values = self.q_values_pileup.to_dense_pileup(
            # self.graph.node_indexes[-1])
        # self.q_values_pileup = q_values

    def call_peaks_from_q_values(self, experiment_info, config=None):
        assert self.q_values_pileup is not None
        caller = CallPeaksFromQvalues(
            self.graph, self.q_values_pileup,
            experiment_info, self.out_file_base_name,
            touched_nodes=self.touched_nodes,
            config=config
            )
        caller.callpeaks()
        self.max_path_peaks = caller.max_paths

    def save_max_path_sequences_to_fasta_file(self, file_name,
                                              sequence_retriever):
        assert self.max_path_peaks is not None
        f = open(self.out_file_base_name + file_name, "w")
        i = 0
        for max_path in self.max_path_peaks:
            seq = sequence_retriever.get_interval_sequence(max_path)
            f.write(">peak" + str(i) + " " +
                    max_path.to_file_line() + "\n" + seq + "\n")
            i += 1
        f.close()
        logging.info("Wrote max path sequences to fasta file: %s" % (
            self.out_file_base_name + file_name))

    @classmethod
    def run_from_intervals(
            cls, graph, sample_intervals,
            control_intervals=None, experiment_info=None,
            out_file_base_name="", has_control=True,
            linear_map=None, configuration=None, stop_after_p_values=False):
        caller = cls(graph, out_file_base_name)
        caller.sample_intervals = sample_intervals
        caller.control_intervals = control_intervals

        caller.run_pre_callpeaks(has_control, experiment_info,
                                 linear_map, configuration)
        caller.get_p_values()
        caller.p_values_pileup.track_size = graph.node_indexes[-1]
        if stop_after_p_values:
            np.save(out_file_base_name + "touched_nodes.npy",
                    np.array(list(caller.touched_nodes), dtype="int"))
            return caller.p_values_pileup.to_sparse_files(
                out_file_base_name + "pvalues")

        caller.get_p_to_q_values_mapping()
        caller.get_q_values()
        caller.call_peaks_from_q_values(experiment_info, configuration)
        return caller


class CallPeaksFromQvalues(object):
    def __init__(self, graph, q_values_pileup,
                 experiment_info,
                 out_file_base_name="",
                 cutoff=0.1, raw_pileup=None, touched_nodes=None,
                 config=None, q_values_max_path=False):
        self.graph = graph
        self.q_values = q_values_pileup
        self.info = experiment_info
        self.out_file_base_name = out_file_base_name
        self.cutoff = cutoff
        self.raw_pileup = raw_pileup
        # self.graph.assert_correct_edge_dicts()
        self.touched_nodes = touched_nodes
        self.save_tmp_results_to_file = False
        self.graph_is_partially_ordered = False
        self.q_values_max_path = q_values_max_path

        if config is not None:
            self.cutoff = config.p_val_cutoff
            self.graph_is_partially_ordered = config.graph_is_partially_ordered
            self.save_tmp_results_to_file = config.save_tmp_results_to_file

        self.info.to_file(self.out_file_base_name + "experiment_info.pickle")
        logging.info("Using p value cutoff %.4f" % self.cutoff)

    def __threshold(self):
        threshold = -np.log10(self.cutoff)
        logging.info("Thresholding peaks on q value %.4f" % threshold)
        self.pre_processed_peaks = self.q_values.threshold_copy(threshold)

        if self.save_tmp_results_to_file:
            self.pre_processed_peaks.to_sparse_files(
                self.out_file_base_name + "peak_areas_pre_postprocess")
            logging.info("Wrote peak pileup before post processing to file.")

    def __postprocess(self):
        logging.info("Filling small Holes")
        if isinstance(self.pre_processed_peaks, DensePileup):
            self.pre_processed_peaks = SparseValues.from_dense_pileup(
                self.pre_processed_peaks.data._values)
        self.pre_processed_peaks = HolesCleaner(
            self.graph,
            self.pre_processed_peaks,
            self.info.read_length,
            self.touched_nodes
        ).run()
        if self.save_tmp_results_to_file:
            self.pre_processed_peaks.to_sparse_files(
                self.out_file_base_name + "peaks_after_hole_cleaning")
            logging.info("Wrote results after holes cleaning to file")
        else:
            logging.info("Not writing results after hole cleaning")

        logging.info("Removing small peaks")

        if isinstance(self.pre_processed_peaks, DensePileup) or isinstance(self.pre_processed_peaks, SparseValues):
            # If dense pileup, we are filtering small peaks while trimming later
            self.filtered_peaks = self.pre_processed_peaks
            logging.info("Not removing small peaks.")
        else:
            self.filtered_peaks = self.pre_processed_peaks.remove_small_peaks(
                self.info.fragment_length)
            logging.info("Small peaks removed")

    def trim_max_path_intervals(self, intervals, end_to_trim=-1, trim_raw=False):
        # Right trim max path intervals, remove right end where q values are 0
        # If last base pair in interval has 0 in p-value, remove hole size
        logging.info("Trimming max path intervals. End: %d" % end_to_trim)
        new_intervals = []
        n_intervals_trimmed = 0
        intervals_too_short = []
        intervals_too_short_untrimmed = []
        for interval in intervals:
            if np.all([rp < 0 for rp in interval.region_paths]):
                use_interval = interval.get_reverse()
                use_interval.score = interval.score
            else:
                use_interval = interval
                assert np.all([rp > 0 for rp in interval.region_paths]), \
                "This method only supports intervals with single rp direction"

            if not trim_raw:
                pileup_values = self.q_values.data.get_interval_values(use_interval)
                pileup_values = 1 * (pileup_values >= -np.log10(self.cutoff))
            else:
                pileup_values = self.raw_pileup.data.get_interval_values(use_interval)
            assert len(pileup_values) == use_interval.length()

            if end_to_trim == 1:
                pileup_values = pileup_values[::-1]

            cumsum = np.cumsum(pileup_values)
            n_zeros_beginning = np.sum(cumsum == 0)

            if n_zeros_beginning == use_interval.length():
                logging.warning("Trimming interval of length %d with %d bp"
                                % (use_interval.length(), n_zeros_beginning))
                continue

            if end_to_trim == -1:
                new_interval = use_interval.get_subinterval(n_zeros_beginning, use_interval.length())
            else:
                new_interval = use_interval.get_subinterval(0, use_interval.length() - n_zeros_beginning)


            if new_interval.length() != use_interval.length():
                n_intervals_trimmed += 1

            new_interval.score = use_interval.score

            if new_interval.length() < self.info.fragment_length:
                intervals_too_short.append(new_interval)
                intervals_too_short_untrimmed.append(use_interval)
                continue
            new_interval = Peak(new_interval.start_position,
                                new_interval.end_position,
                                new_interval.region_paths,
                                new_interval.graph,
                                score=use_interval.score)

            new_intervals.append(new_interval)

        logging.info("N peaks removed because too short after trimming: %d. Written to file." % len(intervals_too_short))
        intervals_too_short = PeakCollection(intervals_too_short)
        intervals_too_short_untrimmed = PeakCollection(intervals_too_short_untrimmed)
        intervals_too_short.to_file(self.out_file_base_name + "too_short_peaks.intervals", text_file=True)
        intervals_too_short_untrimmed.to_file(self.out_file_base_name + "too_short_peaks_untrimmed.intervals",
                                              text_file=True)

        logging.info("Trimmed in total %d intervals" % n_intervals_trimmed)
        logging.info("N intervals left: %d" % len(new_intervals))

        return new_intervals

    def __get_max_paths(self):
        logging.info("Getting maxpaths")
        if not self.q_values_max_path:
            _pileup = SparseValues.from_sparse_files(
                self.out_file_base_name+"direct_pileup")
            self.raw_pileup = _pileup
        else:
            _pileup = SparseValues.from_dense_pileup(
                self.q_values.data._values)
        logging.info("Running Sparse Max Paths")
        max_paths, sub_graphs = SparseMaxPaths(
            self.filtered_peaks, self.graph, _pileup).run()
        PeakCollection(max_paths).to_file(
            self.out_file_base_name + "max_paths_before_length_filtering.intervalcollection",
            text_file=True)
        logging.info("All max paths found")

        # Create dense q
        if not isinstance(self.q_values, DensePileup):
            q_values = DensePileup(self.graph)
            q_values.data._values = self.q_values.to_dense_pileup(
                self.graph.node_indexes[-1])
            self.q_values = q_values

        for max_path in max_paths:
            max_path.set_score(np.max(
                self.q_values.data.get_interval_values(max_path)))
        pairs = list(zip(max_paths, sub_graphs))
        pairs.sort(key=lambda p: p[0].score, reverse=True)
        logging.info("N unfiltered peaks: %s", len(max_paths))
        pairs = [p for p in pairs if
                 p[0].length() >= self.info.fragment_length]
        logging.info("N filtered peaks: %s", len(pairs))
        file_name = self.out_file_base_name + "max_paths.intervalcollection"
        PeakCollection([p[0] for p in pairs]).to_file(
            file_name,
            text_file=True)
        logging.info("Wrote max paths to %s" % file_name)

        np.savez(self.out_file_base_name + "sub_graphs.graphs",
                 **{"peak%s" % i: p[1]._graph for i, p in enumerate(pairs)})
        np.savez(self.out_file_base_name + "sub_graphs.nodeids",
                 **{"peak%s" % i: p[1]._node_ids for i, p in enumerate(pairs)})
        self.max_paths = [p[0] for p in pairs]
        assert max_paths is not None
        # self.max_paths = max_paths

    def __get_subgraphs(self):
        logging.info("Creating subgraphs from peak regions")
        peaks_as_subgraphs = self.filtered_peaks.to_subgraphs()
        logging.info("Writing subgraphs to file")
        peaks_as_subgraphs.to_file(self.out_file_base_name + "peaks.subgraphs")

        logging.info("Found %d subgraphs" % len(peaks_as_subgraphs.subgraphs))
        binary_peaks = peaks_as_subgraphs
        logging.info("Finding max path through subgraphs")
        BCACollection(binary_peaks).to_file(
            self.out_file_base_name + "bcapeaks.subgraphs")
        self.binary_peaks = binary_peaks

    def callpeaks(self):
        logging.info("Calling peaks")
        self.__threshold()
        self.__postprocess()
        # self.__get_subgraphs()
        self.filtered_peaks.to_sparse_files(
            self.out_file_base_name + "peaks_after_postprocessing")
        self.__get_max_paths()

    def save_max_path_sequences_to_fasta_file(self, file_name, sequence_retriever):
        assert self.max_paths is not None, \
                "Max paths has not been found. Run peak calling first."
        assert sequence_retriever is not None
        # assert isinstance(sequence_retriever, vg.sequences.SequenceRetriever)
        f = open(self.out_file_base_name + file_name, "w")
        i = 0
        for max_path in self.max_paths:
            seq = sequence_retriever.get_interval_sequence(max_path)
            f.write(">peak" + str(i) + " " +
                    max_path.to_file_line() + "\n" + seq + "\n")
            i += 1
        f.close()
        logging.info("Wrote max path sequences to fasta file: %s" % (self.out_file_base_name + file_name))

    @staticmethod
    def intervals_to_fasta_file(interval_collection, out_fasta_file_name, sequence_retriever):
        f = open(out_fasta_file_name, "w")
        i = 0
        for max_path in interval_collection.intervals:
            seq = sequence_retriever.get_interval_sequence(max_path)
            f.write(">peak" + str(i) + " " +
                    max_path.to_file_line() + "\n" + seq + "\n")
            i += 1
            if i % 100 == 0:
                logging.info("Writing sequence # %d" % i)
        f.close()

