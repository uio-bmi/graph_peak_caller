import pyvg.util
import pyvg.vg as vg
from graph_peak_caller.callpeaks import CallPeaks
import offsetbasedgraph as obg
data_folder = "graph_peak_caller/dm_test_data/"


def create_graphs():
    pyvg.util.vg_to_offsetbasedgraphs_per_chromosome(
        data_folder + "x.json", data_folder + "obg")


def translate_intervals():
    vg_graph = vg.Graph.create_from_file(
        data_folder+"x.json",
        limit_to_chromosome="chr4")
    ob_graph = vg_graph.to_obg()
    pyvg.util.vg_mapping_file_to_interval_file(
        "intervals_test", vg_graph, data_folder + "sample_reads.json",
        ob_graph)


if __name__ == "__main__":
    translate_intervals()
    exit()
    chromosome = "chr4"
    obg = obg.Graph.from_file(data_folder+"obg%s.tmp" % chromosome)
    interval_file = pyvg.util.vg_mapping_file_to_interval_file(
        "intervals_test", vg_graph, data_folder + "sample_reads.json", ofbg)
    ofbg.to_file("graph.tmp")

    caller = CallPeaks("graph.tmp", interval_file)
    caller.run()

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
