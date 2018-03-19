from pyfaidx import Fasta
from pyvg.sequences import SequenceRetriever
from offsetbasedgraph.graphtraverser import GraphTraverserUsingSequence
from offsetbasedgraph import GraphWithReversals, IntervalCollection
import logging


def find_linear_path_through_chromosome(chromosome, chromend, fasta_file_name, ob_graph_file_name, vg_graph_file_name):
    genome = Fasta(fasta_file_name)
    seq = str(genome[chromosome][0:50818468]).lower()

    logging.info("Creating sequence retriever")
    sequence_retriever = SequenceRetriever.from_vg_json_graph(vg_graph_file_name)

    graph = GraphWithReversals.from_numpy_file(ob_graph_file_name)

    start_nodes = graph.get_first_blocks()
    assert len(start_nodes) == 1, "Found %d start nodes" % start_nodes
    start_node = start_nodes[0]

    traverser = GraphTraverserUsingSequence(graph, seq, sequence_retriever)
    traverser.search_from_node(start_node)
    path = traverser.get_interval_found()
    path = IntervalCollection(path)
    path.to_file("22_path.intervalcollection", text_file=True)
    logging.info("Done")


find_linear_path_through_chromosome("chr22", 50818468, "../hg19.fasta", "tests/whole_genome/22.nobg", "tests/whole_genome/22.json")

