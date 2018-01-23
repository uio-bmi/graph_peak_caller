import offsetbasedgraph as obg
from .experiment_info import ExperimentInfo
from .callpeaks import CallPeaks, Configuration
from .pvalues import PToQValuesMapper
from .densepileup import DensePileup


class MultipleGraphsCallpeaks:

    def __init__(self, graph_names, graph_file_names,
                 sample_inteval_collections,
                 control_interval_collections, linear_maps,
                 experiment_info,
                 out_base_name="", has_control=False):

        self._config = Configuration(
            skip_read_validation=True, save_tmp_results_to_file=False,
            skip_filter_duplicates=True, p_val_cutoff=0.05,
            graph_is_partially_ordered=True)
        self.run()

    @classmethod
    def from_file_base_names(cls):
        # Create lists of names, graphs, sample, control linear maps
        pass

    def run(self):
        self.run_to_p_values()
        self.create_joined_q_value_mapping()
        self.run_from_p_values()

    def run_to_p_values(self):
        for name, samples, controls, lin_map in zip(
                self.names, self.sample_intervals,
                self.control_intervals, self.linear_maps):
            ob_graph = obg.GraphWithReversals(name)
            graph_size = sum(
                block.length() for block in ob_graph.blocks.values())
            info = ExperimentInfo(
                graph_size, self._fragment_length, self.read_length)
            CallPeaks.run_from_intervals(
                ob_graph, samples, controls, info, name, self._has_control,
                lin_map, self._configuration,
                stop_after_p_values=True)

    def create_joined_q_value_mapping(self):
        self._q_mapper = PToQValuesMapper.from_files(self._base_name)

    def run_from_p_values(self):
        for i, name in enumerate(self.names):
            ob_graph = obg.GraphWithReversals(name)
            caller = CallPeaks(ob_graph, name)
            caller.p_to_q_values_mapping = self._q_mapper
            caller.p_values_pileup = DensePileup.from_sparse_files(name)
            caller = CallPeaks.run_from_p_value_and_mapping_files()
            caller.get_p_to_q_values_mapping()
            caller.get_q_values()
