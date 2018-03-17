import logging
import offsetbasedgraph as obg

from .experiment_info import ExperimentInfo
from .control import SparseControl
from .sample import SamplePileupGenerator


class SampleAndControlCreator:
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

        assert isinstance(sample_intervals, obg.IntervalCollection) \
               or isinstance(sample_intervals, str), \
                "Samples intervals must be either interval collection or a file name"

        assert isinstance(control_intervals, obg.IntervalCollection) \
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
                self.graph, self.sample_intervals, self.control_intervals)
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

    def remove_alignments_not_in_graph(self, intervals, is_control=False):
        for interval in self._get_intervals_in_graph(intervals):
            if interval is not False:
                yield interval

    def filter_duplicates_and_count_intervals(self, intervals, is_control=False):
        interval_hashes = {}
        n_duplicates = 0
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
            if next_rp in self.graph.adj_list[rp]:
                new_dir = 1
            elif next_rp in self.graph.reverse_adj_list[rp]:
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

    def _get_intervals_in_graph(self, intervals):
        # Returns only those intervals that exist in graph
        for interval in intervals:
            self._assert_interval_is_valid(interval)
            if interval.region_paths[0] in self.graph.blocks:
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
        else:
            logging.info("Scaling control pileup down with ratio %.3f" % ratio)
            self._control_pileup *= ratio

    def find_info(self):
        sizes = (block.length() for block in self.graph.blocks.values())

        self.genome_size = sum(sizes)
        self.n_reads = sum(1 for line in open(self.control_file_name))

    def create_control(self):
        logging.info("Creating control track using linear map %s" % self.linear_map)
        extensions = [self.info.fragment_length, 1000, 10000] if self.has_control else [10000]
        sparse_control = SparseControl(
            self.linear_map, self.graph, extensions,
            self.info.fragment_length, self.touched_nodes)
        if self.use_global_min_value:
            sparse_control.set_min_value(self.use_global_min_value)
        control_pileup = sparse_control.create(self.control_intervals)
        self.linear_map = None

        # control_pileup.graph = self.ob_graph
        logging.info("Number of control reads: %d" % self.info.n_control_reads)
        self._control_pileup = control_pileup

        # Delete control pileup
        self.control_intervals = None

    def create_sample_pileup(self):
        logging.info("Creating sample pileup")
        logging.info(self.sample_intervals)
        generator = SamplePileupGenerator(
            self.graph, self.info.fragment_length-self.info.read_length)
        pileup = generator.run(self.sample_intervals,
                               self.out_file_base_name+"direct_pileup")
        self.touched_nodes = pileup.touched_nodes

        self._sample_track = self.out_file_base_name + "sample_track"
        self._sample_pileup = pileup
        logging.info("Found in total %d duplicate reads that were removed"
                     % self.n_duplicates_found_sample)

        # Delete sample intervals
        self.sample_intervals = None
