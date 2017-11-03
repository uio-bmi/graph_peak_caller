import pyvg.util
import pyvg.vg as vg
from graph_peak_caller.callpeaks import CallPeaks
import offsetbasedgraph as obg
data_folder = "graph_peak_caller/dm_test_data/"


from graph_peak_caller.util import sparse_maximum, sanitize_indices_and_values
import numpy as np
import sys

#indices, values = sanitize_indices_and_values(np.array([1, 1, 4, 4]), np.array([4, 4, 2, 2]))
#print(indices)
#print(values)
#sys.exit()


indices1 = np.array([0, 4])
values1 = np.array([0, 3])
indices2 = np.array([0, 4])
values2 = np.array([1, 2])


max_indices, max_values = sparse_maximum(indices1, values1, indices2, values2, 10)
print(max_indices)
print(max_values)




def create_graphs():
    pyvg.util.vg_to_offsetbasedgraphs_per_chromosome(
        data_folder + "x.json", data_folder + "obg")


def translate_intervals(interval_file, limimt_to_chromosome=None):
    vg_graph = vg.Graph.create_from_file(
        data_folder+"x.json",
        limit_to_chromosome=limimt_to_chromosome)
    ob_graph = vg_graph.get_offset_based_graph()
    pyvg.util.vg_mapping_file_to_interval_file(
                        data_folder + "intervals_" + limimt_to_chromosome,
                        vg_graph,
                        data_folder + interval_file,
                        ob_graph)


if __name__ == "__main__":
    pass
    #create_graphs()
    #chromosome = "chr2R"
    #translate_intervals("reads3.json", chromosome)

    #obg = obg.Graph.from_file(data_folder+"obg%s.tmp" % chromosome)
    #caller = CallPeaks(data_folder + "obg%s.tmp" % chromosome,
    #                   data_folder + "intervals_" + chromosome, verbose=True)
    #caller.run()

    #caller.create_graph()
    #
    #caller.find_info()
    #caller.determine_shift()
    #caller.sample_file_name = caller.remove_alignments_not_in_graph(caller.sample_file_name)
    #caller.sample_file_name = caller.filter_duplicates(caller.sample_file_name)


#    print("Creating vg graph")
#    # vg_graph = vg.Graph.create_from_file(
#    #"graph_peak_caller/dm_test_data/x.json",
#    #    limit_to_chromosome="chr4",
#    #    do_read_paths=False)
#
#    # vg_graph.to_file("tmp.vggraph")
#    vg_graph = vg.Graph.from_file("tmp.vggraph")
#    print("Creating obg graph")
#    #ob_graph = vg_graph.get_offset_based_graph()
#    # ob_graph.to_file("tmp.obgraph")
#    ob_graph = offsetbasedgraph.Graph.from_file(
#        "tmp.obgraph")
#    print("Reading alignments")
#    f = open("graph_peak_caller/dm_test_data/mapped_reads_sample.json")
#    jsons = (json.loads(line) for line in f.readlines())
#    alignments = [vg.Alignment.from_json(json_dict) for json_dict in jsons]
#    alignments = vg_graph.filter(alignments)
#    obg_alignments = [alignment.path.to_obg(ob_graph)
#                      for alignment in alignments]
#    for alignment in obg_alignments:
#        print(alignment)
#    print("Creating pileup")
#
#
#    pileup = Pileup(ob_graph, obg_alignments, shift=50)
#    pileup.create()
#    pileup.to_bed_graph("tmp.bdg")
#    print("#", pileup.summary())
#     print("#", (sum(interval.length() for interval in obg_alignments)))
