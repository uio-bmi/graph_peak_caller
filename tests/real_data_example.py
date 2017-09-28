from graph_peak_caller import callpeaks
import pyvg
from pyvg.util import vg_gam_file_to_interval_collection
import offsetbasedgraph as obg

def run_with_gam(gam_file_name, vg_graph_file_name, limit_to_chromosomes=False):
    vg_graph = pyvg.Graph.create_from_file(vg_graph_file_name)
    ob_graph = vg_graph.get_offset_based_graph()
    ob_graph.to_file("obgraph")
    #ob_graph = obg.Graph.from_file("obgraph")

    print(ob_graph.adj_list[68566])
    print(ob_graph.adj_list[68567])

    #print(ob_graph.blocks)
    reads_intervals = vg_gam_file_to_interval_collection(None, gam_file_name, ob_graph, max_intervals=100000)
    caller = callpeaks.CallPeaks(ob_graph, reads_intervals, out_file_base_name="real_data_")
    caller.run()


if __name__ == "__main__":
    #run_with_gam("reads3_large.gam", "../graph_peak_caller/dm_test_data/x.json", limit_to_chromosomes=["chr3R", "chr3L", "chr2R", "chrX", "chr2L", "chrY", "chr4"])
    run_with_gam("ENCFF291CUA.gam", "cactus-mhc.json")