from graph_peak_caller.peakcollection import PeakCollection
from graph_peak_caller.util import bed_intervals_to_graph
import offsetbasedgraph as obg
from offsetbasedgraph.graphtraverser import GraphTraverserUsingSequence
from pyvg.sequences import SequenceRetriever
from pybedtools import BedTool

def create_linear_peaks_from_bed(linear_sequence_fasta_file, peaks_bed_file,
                                 obg_graph_file_name, vg_graph_file_name, start_node,
                                 graph_start_offset, graph_end_offset):

    ob_graph = obg.GraphWithReversals.from_file(obg_graph_file_name)
    search_sequence = open(linear_sequence_fasta_file).read()
    sequence_retriever = SequenceRetriever.from_vg_graph(vg_graph_file_name)
    traverser = GraphTraverserUsingSequence(ob_graph, search_sequence, sequence_retriever)
    traverser.search_from_node(start_node)
    linear_path_interval = traverser.get_interval_found()
    print(linear_path_interval.length())
    print(linear_path_interval.region_paths[0])
    print(linear_path_interval.start_position)
    print(linear_path_interval.end_position)

    return
    linear_peaks = PeakCollection.create_from_linear_intervals_in_bed_file(
                        obg_graph_file_name,
                        linear_path_interval,
                        peaks_bed_file,
                        graph_start_offset,
                        graph_end_offset)

    linear_peaks.to_file("linear_peaks")

create_linear_peaks_from_bed("mhc_cleaned2.fa", "../ENCFF155DHA.bed", "cactus-mhc.obg", "cactus-mhc.vg", 225518, 28510119, 33480577)

#linear_peaks = PeakCollection.create_from_linear_intervals_in_bed_file("..ENCFF155DHA.bed")
#graph_peaks = PeakCollection.from_file("real_data_max_paths")

