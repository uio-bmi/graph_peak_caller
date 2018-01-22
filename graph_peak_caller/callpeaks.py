import logging
import numpy as np
import pickle
from offsetbasedgraph import IntervalCollection, DirectedInterval
import pyvg as vg
import offsetbasedgraph
from graph_peak_caller import get_shift_size_on_offset_based_graph
#from .sparsepileupv2 import SparseControlSample
#from .sparsepileupv2 import SparsePileup
from .densepileup import DensePileup, DenseControlSample, QValuesFinder

from .extender import Extender
from .areas import ValuedAreas, BinaryContinousAreas, BCACollection
from .peakscores import ScoredPeak
from .peakcollection import PeakCollection
from . import linearsnarls
IntervalCollection.interval_class = DirectedInterval
from .subgraphcollection import SubgraphCollectionPartiallyOrderedGraph
from .peakcollection import Peak
from memory_profiler import profile
from .pvalues import PValuesFinder, PToQValuesMapper, QValuesFinder
from .experiment_info import ExperimentInfo
from .sampleandcontrolcreator import SampleAndControlCreator

class Configuration:
    def __init__(self, skip_filter_duplicates=False,
                       graph_is_partially_ordered=False,
                       skip_read_validation=False,
                       save_tmp_results_to_file=False,
                       p_val_cutoff=0.1
                ):
        self.skip_filter_duplicates = skip_filter_duplicates
        self.graph_is_partially_ordered = graph_is_partially_ordered
        self.skip_read_validation = skip_read_validation
        self.save_tmp_results_to_file = save_tmp_results_to_file
        self.p_val_cutoff = p_val_cutoff

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

    def run_pre_callpeaks(self, has_control, experiment_info, linear_map, configuration=None):
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
        self.p_values_pileup = \
            PValuesFinder(self.sample_pileup, self.control_pileup).get_p_values_pileup()

    def get_p_to_q_values_mapping(self):
        assert self.p_values_pileup is not None
        finder = PToQValuesMapper.from_p_values_dense_pileup(self.p_values_pileup)
        self.p_to_q_values_mapping = finder.get_p_to_q_values()

    def get_q_values(self):
        assert self.p_values_pileup is not None
        assert self.p_to_q_values_mapping is not None
        finder = QValuesFinder(self.p_values_pileup, self.p_to_q_values_mapping)
        self.q_values_pileup = finder.get_q_values()

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

    @classmethod
    def run_from_intervals(cls, graph, sample_intervals,
                 control_intervals=None, experiment_info=None,
                 out_file_base_name="", has_control=True,
                 linear_map=None, configuration=None):
        caller = cls(graph, out_file_base_name)
        caller.sample_intervals = sample_intervals
        caller.control_intervals = control_intervals

        caller.run_pre_callpeaks(has_control, experiment_info,
                                 linear_map, configuration)
        caller.get_p_values()
        caller.get_p_to_q_values_mapping()
        caller.get_q_values()
        caller.call_peaks_from_q_values(experiment_info, configuration)
        return caller

    def save_max_path_sequences_to_fasta_file(self, file_name, sequence_retriever):
        assert self.max_path_peaks is not None
        f = open(self.out_file_base_name + file_name, "w")
        i = 0
        for max_path in self.max_path_peaks:
            seq = sequence_retriever.get_interval_sequence(max_path)
            f.write(">peak" + str(i) + " " +
                    max_path.to_file_line() + "\n" + seq + "\n")
            i += 1
        f.close()
        logging.info("Wrote max path sequences to fasta file: %s" % (self.out_file_base_name + file_name))


class CallPeaksFromQvalues(object):
    def __init__(self, graph, q_values_pileup,
                 experiment_info,
                 out_file_base_name="",
                 cutoff=0.1, raw_pileup=None, touched_nodes=None,
                 config=None):
        self.graph = graph
        self.q_values = q_values_pileup
        self.info = experiment_info
        self.out_file_base_name = out_file_base_name
        self.cutoff = cutoff
        self.raw_pileup = raw_pileup
        #self.graph.assert_correct_edge_dicts()
        self.touched_nodes = touched_nodes
        self.save_tmp_results_to_file = True
        self.graph_is_partially_ordered = False

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
            self.pre_processed_peaks.to_bed_file(
                self.out_file_base_name + "pre_postprocess.bed")

    def __postprocess(self):
        logging.info("Filling small Holes")

        if isinstance(self.pre_processed_peaks, DensePileup):
            self.pre_processed_peaks.fill_small_wholes_on_dag(
                                    self.info.read_length)
        else:
            self.pre_processed_peaks.fill_small_wholes(
                                    self.info.read_length,
                                    self.out_file_base_name + "_holes.intervals",
                                    touched_nodes=self.touched_nodes)

        logging.info("Removing small peaks")

        # This is slow:
        #self.pre_processed_peaks.to_bed_file(
        #    self.out_file_base_name + "_before_small_peaks_removal.bed")
        #logging.info("Preprocessed peaks written to bed file")


        if isinstance(self.pre_processed_peaks, DensePileup):
            # If dense pileup, we are filtering small peaks while trimming later
            self.filtered_peaks = self.pre_processed_peaks
        else:
            self.filtered_peaks = self.pre_processed_peaks.remove_small_peaks(
                self.info.fragment_length)
        logging.info("Small peaks removed")

    def trim_max_path_intervals(self, intervals, end_to_trim=-1):
        # Right trim max path intervals, remove right end where q values are 0
        # If last base pair in interval has 0 in p-value, remove hole size
        logging.info("Trimming max path intervals. End: %d" % end_to_trim)
        new_intervals = []
        n_intervals_trimmed = 0
        for interval in intervals:
            if np.all([rp < 0 for rp in interval.region_paths]):
                use_interval = interval.get_reverse()
                use_interval.score = interval.score
            else:
                use_interval = interval
                assert np.all([rp > 0 for rp in interval.region_paths]), \
                "This method only supports intervals with single rp direction"

            pileup_values = self.q_values.data.get_interval_values(use_interval)
            assert len(pileup_values) == use_interval.length()

            if end_to_trim == 1:
                pileup_values = pileup_values[::-1]

            cumsum = np.cumsum(pileup_values)
            n_zeros_beginning = np.sum(cumsum == 0)

            if end_to_trim == -1:
                new_interval = use_interval.get_subinterval(n_zeros_beginning, use_interval.length())
            else:
                new_interval = use_interval.get_subinterval(0, use_interval.length() - n_zeros_beginning)

            if new_interval.length() != use_interval.length():
                n_intervals_trimmed += 1
                #print("Trimmed interval into: ")
                #print("   %s" % interval)
                #print("   %s" % new_interval)

            new_interval.score = use_interval.score

            if new_interval.length() < self.info.fragment_length:
                #print("Not keeping too short interval: %s" % new_interval)
                continue
            new_interval = Peak(new_interval.start_position,
                                new_interval.end_position,
                                new_interval.region_paths,
                                new_interval.graph,
                                score=use_interval.score)

            new_intervals.append(new_interval)

        logging.info("Trimmed in total %d intervals" % n_intervals_trimmed)
        logging.info("N intervals left: %d" % len(new_intervals))
        return new_intervals

    def __get_max_paths(self):
        logging.info("Getting maxpaths")
        _pileup = self.raw_pileup if self.raw_pileup is not None else self.q_values
        scored_peaks = (ScoredPeak.from_peak_and_numpy_pileup(peak, _pileup)
                        for peak in self.binary_peaks)
        max_paths = [peak.get_max_path() for peak in scored_peaks]

        max_paths.sort(key=lambda p: p.score, reverse=True)

        if isinstance(self.q_values, DensePileup):
            max_paths = self.trim_max_path_intervals(max_paths, end_to_trim=-1)
            max_paths = self.trim_max_path_intervals(max_paths, end_to_trim=1)


        file_name = self.out_file_base_name + "max_paths.intervalcollection"
        PeakCollection(max_paths).to_file(
            file_name,
            text_file=True)
        logging.info("Wrote max paths to %s" % file_name)

        assert max_paths is not None
        self.max_paths = max_paths

    def __get_subgraphs(self):

        if not self.graph_is_partially_ordered:
            logging.info("Creating subgraphs from peak regions")
            peaks_as_subgraphs = self.filtered_peaks.to_subgraphs()
            logging.info("Writing subgraphs to file")
            peaks_as_subgraphs.to_file(self.out_file_base_name + "peaks.subgraphs")

            logging.info("Found %d subgraphs" % len(peaks_as_subgraphs.subgraphs))
            binary_peaks = [BinaryContinousAreas.from_old_areas(peak) for peak in
                            peaks_as_subgraphs]
            logging.info("Finding max path through subgraphs")
            BCACollection(binary_peaks).to_file(
                self.out_file_base_name + "bcapeaks.subgraphs")
            self.binary_peaks = binary_peaks
        else:
            logging.info("Assuming graph is partially ordered!")
            logging.info("Creating subgraphs from peak regions")
            peaks_as_subgraphs = \
                SubgraphCollectionPartiallyOrderedGraph.create_from_pileup(self.graph, self.filtered_peaks)
            #logging.info("Writing subgraphs to file")
            #peaks_as_subgraphs.to_file(self.out_file_base_name + "peaks.subgraphs")

            logging.info("Found %d subgraphs" % len(peaks_as_subgraphs))
            binary_peaks = peaks_as_subgraphs

            if self.save_tmp_results_to_file:
                logging.info("Writing binary continous peaks to file")
                BCACollection(binary_peaks).to_file(
                    self.out_file_base_name + "bcapeaks.subgraphs")
            self.binary_peaks = binary_peaks

    def callpeaks(self):
        logging.info("Calling peaks")
        self.__threshold()
        self.__postprocess()
        self.__get_subgraphs()
        self.filtered_peaks.to_bed_file(
            self.out_file_base_name + "final_peaks.bed")
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


"""
class CallPeaksFromPValues(object):
    def __init__(self, graph, base_name, p_values_pileup, experiment_info, q_values_mapping):
        self.p_values_pileup = p_values_pileup
        self.q_values_mapping = q_values_mapping
        self.info = experiment_info
        self.base_name = base_name

        self.run()

    def get_q_values(self):
        finder = QValuesFinder(self.p_values_pileup, self.q_values_mapping)
        return finder.get_q_values()

    def run(self):
        q_values = self.get_q_values()
        CallPeaksFromQvalues(self.graph, q_values, self.info, self.base_name)

    @classmethod
    def run_from_files(cls, graph_base_name, base_name,
                         q_values_mapping_file_name):
        graph = offsetbasedgraph.GraphWithReversals.from_numpy_files(graph_base_name)
        experiment_info = ExperimentInfo.from_file(base_name + "experiment_info.pickle")
        with open(q_values_mapping_file_name, "rb") as f:
            q_values_mapping = pickle.load(f)
        p_values_pileup = DensePileup.from_sparse_files(graph, base_name)
        return cls(graph, p_values_pileup, experiment_info, q_values_mapping,
                   base_name)
"""


if __name__ == "__main__":
    chromosome = "chr2R"
    vg_graph = vg.Graph.create_from_file(
        "dm_test_data/x_%s.json" % chromosome, 30000, chromosome)
    ofbg = vg_graph.get_offset_based_graph()
    interval_file = vg.util.vg_mapping_file_to_interval_file(
        "intervals_reads3_chr2R", vg_graph,
        "dm_test_data/reads3_small.json", ofbg)
    ofbg.to_file("graph.tmp")

    caller = CallPeaks("graph.tmp", interval_file)
    caller.create_graph()
    caller.find_info()
    caller.determine_shift()
    caller.sample_file_name = caller.remove_alignments_not_in_graph(
        caller.sample_file_name)
    caller.sample_file_name = caller.filter_duplicates(caller.sample_file_name)
