import logging
import pickle
import numpy as np
from offsetbasedgraph import IntervalCollection, DirectedInterval
import offsetbasedgraph
from .densepileup import DensePileup

from .extender import Extender
from .areas import ValuedAreas
from .sparsepvalues import PValuesFinder, PToQValuesMapper
from .experiment_info import ExperimentInfo
# from .directsamplepileup import main as samplemain
from .sparsesampleandcontrolcreator import SparseControl
from .sparsediffs import SparseDiffs
from .sparsegraphpileup import SamplePileupGenerator


IntervalCollection.interval_class = DirectedInterval


def enable_filewrite(func):
    def wrapper(*args, **kwargs):
        intervals = args[1]
        if isinstance(intervals, str):
            intervals = IntervalCollection.from_file(intervals)

        write_to_file = kwargs.pop("write_to_file", False)
        interval_list = func(args[0], intervals, **kwargs)

        if write_to_file:
            interval_collection = IntervalCollection(interval_list)
            interval_collection.to_file(write_to_file)
            return write_to_file
        else:
            return interval_list

    return wrapper


class SampleAndControlCreatorO(object):
    def __init__(self, graph, sample_intervals,
                 control_intervals=None, experiment_info=None,
                 verbose=False, out_file_base_name="", has_control=True,
                 linear_map=None, skip_filter_duplicates=False,
                 skip_read_validation=False,
                 save_tmp_results_to_file=True,
                 configuration=None,
                 graph_is_partially_ordered=False):
        """
        :param sample_intervals: Either an interval collection or file namen
        :param control_intervals: Either an interval collection or a file name
        """

        assert linear_map is not None, "LinearMap cannot be None"
        assert isinstance(linear_map, str), "Must be file name"

        assert isinstance(sample_intervals, IntervalCollection) \
               or isinstance(sample_intervals, str), \
                "Samples intervals must be either interval collection or a file name"

        assert isinstance(control_intervals, IntervalCollection) \
               or isinstance(control_intervals, str), \
                "control_intervals must be either interval collection or a file name"

        if not has_control:
            logging.info("Running without control")
        else:
            logging.info("Running with control")

        self.graph = graph

        self.sample_intervals = sample_intervals
        self.control_intervals = control_intervals
        self.has_control = has_control
        self.linear_map = linear_map

        self._p_value_track = "p_value_track"
        self._q_value_track = "q_value_track"
        self.info = experiment_info
        self.verbose = verbose
        self._control_pileup = None
        self._sample_pileup = None
        self.out_file_base_name = out_file_base_name
        self.pre_processed_peaks = None
        self.filtered_peaks = None
        self.skip_filter_duplicates = skip_filter_duplicates
        self.skip_read_validation = skip_read_validation
        self.save_tmp_results_to_file = save_tmp_results_to_file
        self.n_duplicates_found_sample = 0

        self.max_paths = None
        self.peaks_as_subgraphs = None

        self.create_graph()
        self.touched_nodes = None  # Nodes touched by sample pileup
        self.graph_is_partially_ordered = graph_is_partially_ordered

        # Tmp hack, shuold be default
        if configuration is not None:
            self.graph_is_partially_ordered = configuration.graph_is_partially_ordered
            self.skip_filter_duplicates = configuration.skip_filter_duplicates
            self.skip_read_validation = configuration.skip_read_validation
            self.save_tmp_results_to_file = configuration.save_tmp_results_to_file
            self.use_global_min_value = configuration.use_global_min_value


            if self.use_global_min_value is not None:
                logging.info("Will use global min value %.4f when creating background signal" % self.use_global_min_value)

        if self.skip_filter_duplicates:
            logging.warning("Not removing duplicates")
        else:
            logging.info("Will remove duplicates.")

    def run(self):
        self.preprocess()
        if self.info is None:
            self.info = ExperimentInfo.find_info(
                self.ob_graph, self.sample_intervals, self.control_intervals)
        self.create_sample_pileup()
        self.create_control()
        self.scale_tracks()

    def preprocess(self):
        self.info.n_control_reads = 0
        self.info.n_sample_reads = 0

        if not self.skip_read_validation:
            self.sample_intervals = self.remove_alignments_not_in_graph(
                                        self.sample_intervals)
        else:
            logging.warning("Skipping validation of reads. Not checking whether reads are valid"
                            " or inside the graph.")

        self.sample_intervals = self.filter_duplicates_and_count_intervals(
                                    self.sample_intervals, is_control=False)

        if not self.skip_read_validation:
            self.control_intervals = self.remove_alignments_not_in_graph(
                                    self.control_intervals, is_control=True)
        self.control_intervals = self.filter_duplicates_and_count_intervals(
                                    self.control_intervals, is_control=True)

    @enable_filewrite
    def remove_alignments_not_in_graph(self, intervals, is_control=False):
        for interval in self._get_intervals_in_ob_graph(intervals):
            if interval is not False:
                yield interval

    @enable_filewrite
    def filter_duplicates_and_count_intervals(self, intervals, is_control=False):
        interval_hashes = {}
        n_duplicates = 0
        n_reads_left = 0
        for interval in intervals:
            if not self.skip_filter_duplicates:
                hash = interval.hash(ignore_end_pos=True)
                if hash in interval_hashes:
                    n_duplicates += 1
                    if not is_control:
                        self.n_duplicates_found_sample += 1
                    continue

                interval_hashes[hash] = True

            if is_control:
                self.info.n_control_reads += 1
            else:
                self.info.n_sample_reads += 1

            yield interval

    def _assert_interval_is_valid(self, interval):
        # Checks that the interval (read) is a valid connected interval
        direction = None
        for i, rp in enumerate(interval.region_paths[:-1]):
            next_rp = interval.region_paths[i+1]
            if next_rp in self.ob_graph.adj_list[rp]:
                new_dir = 1
            elif next_rp in self.ob_graph.reverse_adj_list[rp]:
                new_dir = -1
            else:
                logging.error("Invalid interval: Rp %d of interval %s is not "
                              "connected in graph to rp %d, which is the next rp"
                              % (rp, interval, next_rp))
                raise Exception("Invalid interval")

            if direction is None:
                direction = new_dir
            else:
                if new_dir != direction:
                    logging.error("Invalid read: Interval %s is going edges in multiple directions.")
                    raise Exception("Invalid interval")

        return True

    def _get_intervals_in_ob_graph(self, intervals):
        # Returns only those intervals that exist in graph
        for interval in intervals:
            self._assert_interval_is_valid(interval)
            if interval.region_paths[0] in self.ob_graph.blocks:
                yield interval
            else:
                logging.warning("Interval: %s" % interval)
                raise Exception("Interval not in graph")

    def scale_tracks(self, update_saved_files=False):
        logging.info("Scaling tracks to ratio: %d / %d" % (
            self.info.n_sample_reads,
            self.info.n_control_reads))
        ratio = self.info.n_sample_reads/self.info.n_control_reads

        if self.info.n_sample_reads == self.info.n_control_reads:
            logging.info("Not scaling any tracks because of same amount of reads")
            self._control_pileup.to_bed_graph(
                self.out_file_base_name + "scaled_control.bdg")
            return

        if ratio > 1:
            logging.warning("More reads in sample than in control")
            self._sample_pileup *= ratio
            self._sample_pileup.scale(1/ratio)
            self._sample_pileup.to_bed_graph(
                self.out_file_base_name + "scaled_treat.bdg")
            if update_saved_files:
                self._sample_pileup.to_bed_graph(self._sample_track)
        else:
            logging.info("Scaling control pileup down using ration %.3f" % ratio)
            self._control_pileup *= ratio  # .scale(ratio)
            self._control_pileup.to_bed_graph(
                self.out_file_base_name + "scaled_control.bdg")
            if update_saved_files:
                self._control_pileup.to_bed_graph(self._control_track)

    def find_info(self):
        sizes = (block.length() for block in self.ob_graph.blocks.values())

        self.genome_size = sum(sizes)
        self.n_reads = sum(1 for line in open(self.control_file_name))

    def create_graph(self):
        logging.info("Creating graph")
        if isinstance(self.graph, str):
            self.ob_graph = offsetbasedgraph.Graph.from_file(self.graph)
        else:
            self.ob_graph = self.graph
            logging.info("Graph already created")

    def create_control(self):
        logging.info("Creating control track using linear map %s" % self.linear_map)
        extensions = [self.info.fragment_length, 1000, 10000] if self.has_control else [10000]
        sparse_control = SparseControl(
            self.linear_map, self.graph, extensions,
            self.info.fragment_length, self.touched_nodes)
        if self.use_global_min_value:
            sparse_control.set_min_value(self.use_global_min_value)
        control_pileup = sparse_control.create(self.control_intervals)
        control_pileup.to_sparse_files(self.out_file_base_name + "control_track")
        # control_pileup = SparseDiffs.from_dense_pileup(
        #     control_pileup.data._values)
        # control_pileup = linearsnarls.create_control(
        #     self.linear_map,  self.control_intervals,
        #     extensions, self.info.fragment_length,
        #     ob_graph=self.graph,
        #     touched_nodes=self.touched_nodes,
        #     use_global_min_value=self.use_global_min_value
        # )

        # Delete linear map
        self.linear_map = None

        # control_pileup.graph = self.ob_graph
        logging.info("Number of control reads: %d" % self.info.n_control_reads)

        if self.save_tmp_results_to_file and False:
            self._control_track = self.out_file_base_name + "control_track.bdg"
            control_pileup.to_bed_graph(self._control_track)
            logging.info("Saved control pileup to " + self._control_track)

        self._control_pileup = control_pileup

        # Delete control pileup
        self.control_intervals = None

    def get_score(self):
        logging.info("Getting p valyes.")
        p_values_finder = PValuesFinder(self._sample_pileup, self._control_pileup)
        self._p_values_pileup = p_values_finder.get_p_values_pileup()

        # Delete sample and control pileups
        self._control_pileup = None
        self._sample_pileup = None

    def get_p_to_q_values_mapping(self):
        return PToQValuesMapper.from_p_values_dense_pileup(self.p_values)

    #@profile
    def create_sample_pileup(self):
        logging.debug("In sample pileup")
        logging.info("Creating sample pileup")
        alignments = self.sample_intervals
        logging.info(self.sample_intervals)
        extender = Extender(self.ob_graph, self.info.fragment_length)
        valued_areas = ValuedAreas(self.ob_graph)
        logging.info("Extending sample reads")
        areas_list = (extender.extend_interval(interval)
                      for interval in alignments)
        i = 0
        logging.info("Processing areas")

        #touched_nodes = set()  # Speedup thing, keep track of nodes where areas are on
        pileup = DensePileup.create_from_binary_continous_areas(
                    self.ob_graph, areas_list)
        # print(np.sum(pileup.data._values))
        # print(len(pileup.data._touched_nodes))
        for node_id in pileup.data._touched_nodes:
            assert np.sum(pileup.data.values(node_id)) > 0

        touched_nodes = pileup.data._touched_nodes
        self.touched_nodes = touched_nodes

        self._sample_track = self.out_file_base_name + "sample_track.bdg"
        if self.save_tmp_results_to_file:
            logging.info("Saving sample pileup to file")
            pileup.to_bed_graph(self._sample_track)
            logging.info("Saved sample pileup to " + self._sample_track)

            logging.info("Writing touched nodes to file")
            with open(self.out_file_base_name + "touched_nodes.pickle", "wb") as f:
                pickle.dump(touched_nodes, f)

            logging.info("N touched nodes: %d" % len(touched_nodes))
        else:
            logging.info("Not saving sample pileup to files.")

        self._sample_pileup = pileup

        # Delete sample intervals
        self.sample_intervals = None

    def _write_vg_alignments_as_intervals_to_bed_file(self):

        pass


class SampleAndControlCreator(SampleAndControlCreatorO):
    def create_sample_pileup(self):
        logging.debug("In sample pileup")
        logging.info("Creating sample pileup")
        logging.info(self.sample_intervals)
        logging.info("Processing areas")
        generator = SamplePileupGenerator(
            self.graph, self.info.fragment_length-self.info.read_length)
        print("#######", self.out_file_base_name+"direct_pileup")
        pileup = generator.run(self.sample_intervals,
                               self.out_file_base_name+"direct_pileup")
        # pileup = samplemain(self.sample_intervals, self.graph,
        #                     self.info.fragment_length-self.info.read_length,
        #                     self.out_file_base_name+"direct_pileup")
        # assert np.all(pileup.data._values == pileup2.data._values)
        self.touched_nodes = pileup.touched_nodes

        self._sample_track = self.out_file_base_name + "sample_track"
        if self.save_tmp_results_to_file or True:
            logging.info("Saving sample pileup to file")
            pileup.to_sparse_files(self._sample_track)
            logging.info("Saved sample pileup to " + self._sample_track)

            logging.info("Writing touched nodes to file")
            with open(self.out_file_base_name + "touched_nodes.pickle", "wb") as f:
                pickle.dump(self.touched_nodes, f)

            logging.info("N touched nodes: %d" % len(self.touched_nodes))
        else:
            logging.info("Not saving sample pileup to files.")

        self._sample_pileup = pileup
        logging.info("Found in total %d duplicate reads that were removed" % self.n_duplicates_found_sample)

        # Delete sample intervals
        self.sample_intervals = None
