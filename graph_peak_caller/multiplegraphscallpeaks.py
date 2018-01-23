from .callpeaks import  CallPeaks


class MultipleGraphsCallpeaks:

    def __init__(self, graph_names, graph_file_names,
                 sample_file_names,
                 control_file_names, linear_maps,
                 experiment_info,
                 out_base_name="", has_control=False):

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
        for i, name in enumerate(self.names):
            caller = CallPeaks.run_from_intervals(stop_after_p_values=True)

    def create_joined_q_value_mapping(self):
        pass

    def run_from_p_values(self):
        for i, name in enumerate(self.names):
            caller = CallPeaks.run_from_p_value_and_mapping_files()


