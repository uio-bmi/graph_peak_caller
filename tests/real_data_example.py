from graph_peak_caller import callpeaks
import pyvg
from pyvg.util import vg_gam_file_to_intervals


def run_with_gam_(gam_file_name, vg_graph_file_name):
    vg_graph = pyvg.Graph.from_file(vg_graph_file_name)
    ob_graph = vg_graph.get_offset_based_graph()
    reads_intervals = vg_gam_file_to_intervals(vg_graph, gam_file_name, ob_graph)
    caller = callpeaks.CallPeaks(ob_graph, reads_intervals, out_file_base_name="real_data_")
    caller.run()


if __name__ == "__main__":
    run_with_gam("../graph_peak_caller/dm_test_data/....gam", "../graph_peak_caller/dm_test_data/graph.vg")