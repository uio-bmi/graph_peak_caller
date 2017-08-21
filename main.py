import graph_peak_caller.vg as vg
import offsetbasedgraph
import json
from graph_peak_caller.pileup import Pileup


if __name__ == "__main__":
    print("Creating vg graph")
    # vg_graph = vg.Graph.create_from_file(
    #"graph_peak_caller/dm_test_data/x.json",
    #    limit_to_chromosome="chr4",
    #    do_read_paths=False)

    # vg_graph.to_file("tmp.vggraph")
    vg_graph = vg.Graph.from_file("tmp.vggraph")
    print("Creating obg graph")
    #ob_graph = vg_graph.get_offset_based_graph()
    # ob_graph.to_file("tmp.obgraph")
    ob_graph = offsetbasedgraph.Graph.from_file(
        "tmp.obgraph")
    print("Reading alignments")
    f = open("graph_peak_caller/dm_test_data/mapped_reads_sample.json")
    jsons = (json.loads(line) for line in f.readlines())
    alignments = [vg.Alignment.from_json(json_dict) for json_dict in jsons]
    alignments = vg_graph.filter(alignments)
    obg_alignments = [alignment.path.to_obg(ob_graph)
                      for alignment in alignments]
    for alignment in obg_alignments:
        print(alignment)
    print("Creating pileup")


    pileup = Pileup(ob_graph, obg_alignments, shift=50)
    pileup.create()
    pileup.to_bed_graph("tmp.bdg")
    print("#", pileup.summary())
    print("#", (sum(interval.length() for interval in obg_alignments)))
