from collections import defaultdict
import json
from offsetbasedgraph import graph, Interval, Translation
from Bio import SeqIO
import re
import pickle

def create_graph_from_vg_json(file_name):
    """
    Creates an offsetbasedgraph and a dict of block
    sequences (block id: sequences)
    """
    f = open(file_name)
    block_sequences = {}

    lines = f.readlines()
    n_lines = len(lines)
    print("Number of lines: %d" % n_lines)
    blocks_offset_based_graph = {}
    edges_offset_based_graph = defaultdict(list)
    paths = defaultdict(list) # path name (chromosome) : mapping


    i = 0
    for line in lines:
        print("Line: %d/%d %s" % (i , n_lines, line[0:100]))
        #print("Line: %d/%d" % (i , n_lines))

        line = json.loads(line)

        if "path" in line:
            p = line["path"][0]
            paths[p["name"]].append(p["mapping"])

            print("Path %s" % p["name"])
            #print("Mapping: \n%s" % p["mapping"])
            #print(line["node"])


        if i >= 100:
            break

        i += 1

        #    print(line["path"])


        if not "node" in line:
            print("No nodes in line")
            continue

        nodes = line["node"]
        # print("Number of nodes: %d" % len(nodes))

        for node in nodes:
            id = node["id"]
            length = len(node["sequence"])
            # print(id, length)
            blocks_offset_based_graph[id] = length
            block_sequences[id] = node["sequence"].lower()

        if "edge" in line:
            edges = line["edge"]

            for edge in edges:
                edges_offset_based_graph[edge["from"]].append(edge["to"])
        else:
            print("Found no edges")

    g = graph.Graph(blocks_offset_based_graph, edges_offset_based_graph)

    return g, block_sequences


def create_translation_from_vg_json(json_file_name, g, block_sequences):
    chromosome, chromosome_sequence = map_subgraph_to_chromosome(
        "dm6.fa", g, block_sequences)
    path = find_block_through_graph_matching_sequence(
        chromosome_sequence, g, block_sequences)

    linear_blocks = {}
    linear_blocks[chromosome] = len(chromosome_sequence)

    linear_graph = graph.Graph(linear_blocks, {})
    last_block = path[-1]
    last_block_length = len(block_sequences[last_block])
    path_interval = Interval(0, last_block_length, path)
    forward_trans = {}
    forward_trans[chromosome] = [path_interval]

    backward_trans = {}
    offset = 0
    for block in path:
        block_length = len(block_sequences[block])
        backward_trans[block] = [
            Interval(0, block_length, [chromosome], linear_graph)]
        offset += block_length

    trans = Translation(forward_trans, backward_trans, linear_graph)

    return trans


def map_subgraph_to_chromosome(fasta_file, sub_graph, block_sequences):
    """
    Searches the fasta file for the chromosome
    with starting sequence that matches the beginning of the graph.
    Returns the chromosome name.
    """
    start_blocks = g.get_first_blocks()
    assert len(start_blocks) == 1, "len start blocks != 1"
    start_block = start_blocks[0]
    start_sequence = block_sequences[start_block]

    return find_chromosome_with_sequence(fasta_file, start_sequence)


def find_chromosome_with_sequence(fasta_file, sequence):

    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    for fasta in fasta_sequences:
        name, fasta_sequence = fasta.id, str(fasta.seq)
        print(name)
        if re.match(sequence, fasta_sequence, re.I):  # Case insensitive startswith check
            print("found")
            print(name)
            return name, fasta_sequence

    assert False, "Could not find start sequence %s in fasta file" % sequence
    return False


def find_block_through_graph_matching_sequence(sequence, graph, block_sequences):
    """
    Returns a list of blocks (that are connected) matching the sequence
    """
    sequence = sequence.lower()
    first_block = graph.get_first_blocks()[0]
    path = [first_block]
    current_block = first_block
    offset = len(block_sequences[first_block])
    assert sequence.startswith(block_sequences[first_block]), "%s does not start with %s" % (sequence[0:200], block_sequences[first_block])

    while True:
        edges = graph.adj_list[current_block]
        if len(edges) == 0:
            break

        next_block = None
        for potential_next in edges:
            if potential_next not in graph.blocks:
                return path

            if sequence[offset:].startswith(block_sequences[potential_next]):
                next_block = potential_next
                #print("  Match with %s "  % potential_next)
                break
            #else:
            #    print("  Potential edge %s did not match" % potential_next)

        assert next_block is not None
        #print("Next block: %s" % next_block)

        offset += len(block_sequences[next_block])
        path.append(next_block)
        current_block = next_block

    return path


def save_graph_and_block_sequences(file_name, g, block_sequences):
    g.to_file(file_name + ".offsetbasedgraph")
    pickle_file = open(file_name + ".block_sequences", 'wb')
    pickle.dump(block_sequences, pickle_file)
    pickle_file.close()


def get_graph_and_block_sequences_from_files(file_name):
    g = graph.Graph.from_file(file_name + ".offsetbasedgraph")
    pickle_file = open(file_name + ".block_sequences", 'rb')
    block_sequences = pickle.load(pickle_file)
    pickle_file.close()

    return g, block_sequences


g, block_sequences = create_graph_from_vg_json("dm_test_data/x.json")
#save_graph_and_block_sequences("dm_test_data/x", g, block_sequences)
#g2, block_sequences2 = get_graph_and_block_sequences_from_files("x")

#assert g2 == g, "Not equal"

# trans = create_translation_from_vg_json("x.json", g, block_sequences)

# print(trans)




