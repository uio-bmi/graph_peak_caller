import logging

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

        assert linear_map is not None, "LinearMap cannot be None"
        assert isinstance(linear_map, str), "Must be file name"

        if not has_control:
            logging.info("Running without control")
        else:
            logging.info("Running with control")

        self.graph = graph

        self.has_control = has_control
        self.linear_map = linear_map

        self.info = experiment_info
        self._control_pileup = None
        self._sample_pileup = None
        self.out_file_base_name = out_file_base_name
        self.skip_filter_duplicates = skip_filter_duplicates
        self.save_tmp_results_to_file = save_tmp_results_to_file
        self.n_duplicates_found_sample = 0
        self.touched_nodes = None  # Nodes touched by sample pileup

        # Tmp hack, shuold be default
        if configuration is not None:
            self.skip_filter_duplicates = configuration.skip_filter_duplicates
            self.skip_read_validation = configuration.skip_read_validation
            self.save_tmp_results_to_file = configuration.save_tmp_results_to_file
            self.use_global_min_value = configuration.use_global_min_value

            if self.use_global_min_value is not None:
                logging.info("Will use global min value %.4f when creating background signal" % self.use_global_min_value)
        if self.skip_filter_duplicates:
            self.sample_intervals = Intervals(sample_intervals)
            self.control_intervals = Intervals(control_intervals)
        else:
            self.sample_intervals = UniqueIntervals(sample_intervals)
            self.control_intervals = UniqueIntervals(control_intervals)


    def run(self):
        if self.info is None:
            self.info = ExperimentInfo.find_info(
                self.graph, self.sample_intervals, self.control_intervals)
        self.create_sample_pileup()
        self.create_control()
        self.scale_tracks()

    def scale_tracks(self, update_saved_files=False):
        logging.info("Scaling tracks to ratio: %d / %d" % (
            self.sample_intervals.n_reads,
            self.control_intervals.n_reads))
        ratio = self.sample_intervals.n_reads/self.control_intervals.n_reads
        if ratio == 1:
            print("NO RATIO")
            return

        if ratio > 1:
            logging.warning("More reads in sample than in control")
            self._sample_pileup *= 1/ratio
            # self._sample_pileup.scale(1/ratio)
        else:
            logging.info("Scaling control pileup down with ratio %.3f" % ratio)
            self._control_pileup *= ratio

    def find_info(self):
        sizes = (block.length() for block in self.graph.blocks.values())
        self.genome_size = sum(sizes)
        self.n_reads = sum(1 for line in open(self.control_file_name))

    def create_control(self):
        logging.info("Creating control track using linear map %s"
                     % self.linear_map)
        extensions = [self.info.fragment_length, 1000, 10000] if self.has_control else [10000]
        sparse_control = SparseControl(
            self.linear_map, self.graph, extensions,
            self.info.fragment_length, self.touched_nodes)
        if self.use_global_min_value:
            sparse_control.set_min_value(self.use_global_min_value)
        control_pileup = sparse_control.create(self.control_intervals)
        self.linear_map = None
        self.info.n_control_reads = self.control_intervals.n_reads
        # control_pileup.graph = self.ob_graph
        logging.info(
            "Number of control reads: %d" % self.control_intervals.n_reads)
        self._control_pileup = control_pileup
        print(control_pileup)

    def create_sample_pileup(self):
        logging.info("Creating sample pileup")
        print("---------CREATING SAMPLE_____________")
        generator = SamplePileupGenerator(
            self.graph, self.info.fragment_length-self.info.read_length)
        pileup = generator.run(self.sample_intervals,
                               self.out_file_base_name+"direct_pileup")
        self.touched_nodes = pileup.touched_nodes
        
        self._sample_track = self.out_file_base_name + "sample_track"
        self._sample_pileup = pileup
        print(pileup.get_sparse_values())
        logging.info("N Sample Reads %s" % self.sample_intervals.n_reads)
        logging.info("Found in total %d duplicate reads that were removed"
                     % self.sample_intervals.n_duplicates)
        self.info.n_sample_reads = self.sample_intervals.n_reads


class Intervals:
    def __init__(self, intervals):
        self._intervals = intervals
        self.n_reads = 0
        self.n_duplicates = 0

    def _count_and_return(self, x):
        self.n_reads += 1
        return x

    def __iter__(self):
        for interval in self._intervals:
            self.n_reads += 1
            yield interval


class UniqueIntervals:
    def __init__(self, intervals):
        self._intervals = intervals
        self.n_reads = 0
        self.n_duplicates = 0

    def __iter__(self):
        interval_hashes = {}
        for interval in self._intervals:
            _hash = interval.hash(ignore_end_pos=True)
            if _hash in interval_hashes:
                self.n_duplicates += 1
                continue

            interval_hashes[_hash] = True
            self.n_reads += 1
            yield interval
