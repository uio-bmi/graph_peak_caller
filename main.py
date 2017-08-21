import graph_peak_caller.vg as vg
import offsetbasedgraph
import json
from graph_peak_caller.pileup import Pileup
from graph_peak_caller.get_lambda_track import BackgroundTrack, ControlInfo


class Pipeline(object):

    def __init__(self, graph_file_name, sample_file_name, control_file_name=None, chromosome=None):
        self.graph_file_name = graph_file_name
        self.sample_file_name = sample_file_name
        self.control_file_name = control_file_name if control_file_name is not None else sample_file_name
        self.has_control = control_file_name is not None
        self.create_graphs()
        self.chromosome = chromosome

    def 

    def create_graphs(self):
        self.vg_graph = vg.Graph.create_from_file(self.graph_file_name)
        self.ob_graph = vg_graph.get_offset_based_graph()

    def create_control(self):
        bg_track = BackgroundTrack(self.ob_graph, control_file_name, self.d, 1000, 10000, 
        f = open()

        jsons = (json.loads(line) for line in f.readlines())
        alignments = [vg.Alignment.from_json(json_dict) for json_dict in jsons]        

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
