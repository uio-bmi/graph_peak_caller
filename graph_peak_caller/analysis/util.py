import logging
import offsetbasedgraph as obg
from offsetbasedgraph.graphtraverser import GraphTraverserUsingSequence
import pyvg


class LinearRegion(object):
    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.start = start
        self.end = end

LRC_REGION = LinearRegion("chr19", 54025634, 55084318)
MHC_REGION = LinearRegion("chr6", 28510119, 33480577)


def bed_intervals_to_graph(obg_graph, linear_path_interval,
                           bed_file_name, graph_start_offset):
    f = open(bed_file_name)
    intervals_on_graph = []
    for line in f:
        l = line.split()
        start = int(l[1])
        end = int(l[2])

        start = start - graph_start_offset
        end = end - graph_start_offset
        intervals_on_graph.append(linear_path_interval.get_subinterval(start, end))

    return intervals_on_graph


def fasta_sequence_to_linear_path_through_graph(
        linear_sequence_fasta_file, sequence_retriever, ob_graph, start_node):
    search_sequence = open(linear_sequence_fasta_file).read()
    print("Length of search sequence: %d" % len(search_sequence))
    traverser = GraphTraverserUsingSequence(
        ob_graph, search_sequence, sequence_retriever)
    traverser.search_from_node(start_node)
    linear_path_interval = traverser.get_interval_found()
    return linear_path_interval


def get_linear_paths_in_graph(ob_graph, vg_graph, write_to_file_name=None):
    assert ob_graph is not None
    intervals = {}
    for path in vg_graph.paths:
        obg_interval = path.to_obg(ob_graph=ob_graph)
        obg_interval.name = path.name
        print("Path name: %s" % path.name)
        intervals[obg_interval.name] = obg_interval

    if write_to_file_name is not None:
        logging.info("Writing linear path to %s" % write_to_file_name)
        collection = obg.IntervalCollection(intervals.values())
        collection.to_file(write_to_file_name, text_file=True)

    return intervals


def create_ob_graph_from_vg(vg_json_graph_file_name, ob_graph_file_name="graph.obg"):
    vg_graph = pyvg.Graph.create_from_file(vg_json_graph_file_name)
    ob_graph = vg_graph.get_offset_based_graph()
    ob_graph.to_file(ob_graph_file_name)
    logging.info("Wrote obgraph to %s" % ob_graph_file_name)


def create_linear_path(ob_graph, vg_graph, path_name="ref", write_to_file="linear_path.intervalcollection"):
    assert ob_graph is not None
    linear_paths = get_linear_paths_in_graph(ob_graph, vg_graph, write_to_file_name=write_to_file)
    ref_path = obg.NumpyIndexedInterval.from_interval(linear_paths[path_name])
    return ref_path
