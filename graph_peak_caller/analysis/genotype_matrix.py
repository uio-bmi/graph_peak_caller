import offsetbasedgraph as obg
import numpy as np


class GenoType:
    def __init__(self, start_pos, end_pos, variant_ids):
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.variant_ids = variant_ids


class VariantList:
    def __init__(self, indexed_interval, graph, seq_graph):
        self._indexed_interval = indexed_interval
        self._graph = graph
        self._seq_graph = seq_graph
        self.alts = {}
        self.insertions = {}

    def add_simple_var(self, pos, alt_seq, ref_seq, i):
        for j, pair in enumerate(zip(alt_seq, ref_seq)):
            if pair[0] != pair[1]:
                break
        print(ref_seq, alt_seq)
        alt_seq = alt_seq[j:]
        ref_seq = ref_seq[j:]
        if not alt_seq:
            alt_seq = "."
            
        print(ref_seq, alt_seq)
        pos += j
        assert self._indexed_interval.get_node_offset_at_offset(pos) == 0, (self._indexed_interval.get_node_offset_at_offset(pos), i, j)
        prev_node = self._indexed_interval.get_node_at_offset(pos-1)
        if alt_seq == ".":
            self.insertions[(prev_node, self._indexed_interval.get_node_at_offset(pos+len(ref_seq)))] = i
        else:
            next_nodes = [node for node in self._graph.adj_list[prev_node]
                          if self._seq_graph.get_sequence_on_directed_node(node).lower().startswith(alt_seq.lower())]
            assert len(next_nodes) == 1, (prev_node, next_nodes, i)
            if not next_nodes:
                # print(seq)
                # print([self._seq_graph.get_sequence_on_directed_node(node).lower() for node in self._graph.adj_list[prev_node]])
                return
            next_node = next_nodes[0]
            self.alts[next_node] = i

    def add_var(self, pos, alt_seq, ref_seq, i):
        for seq in alt_seq.split(","):
            self.add_simple_var(pos, seq, ref_seq, i)

    def create_lookup(self, vcf_file_name):
        for i, line in enumerate(open(vcf_file_name)):
            parts = line.split("\t")
            pos = int(parts[1])-1  # to 0-indx
            ref_seq = parts[3]
            alt_seq = parts[4]
            self.add_var(pos, alt_seq, ref_seq, i)

    def create_type(self, interval):
        node_ids = interval.region_paths
        alts = [self.alts[node_id] for node_id in node_ids if node_id in self.alts]
        insertions = [self.insertions[(n1, n2)] for n1, n2 in zip(node_ids[:-1], node_ids[1:])
                      if (n1, n2) in self.insertions]
        start_pos = self._indexed_interval.get_offset_at_position(interval.start_position)
        end_pos = self._indexed_interval.get_offset_at_position(interval.end_position)
        return GenoType(start_pos, end_pos, alts+insertions)

    @classmethod
    def from_name(cls, name):
        indexed_interval = obg.NumpyIndexedInterval.from_file(
            name+"_linear_pathv2.interval")
        graph = obg.Graph.from_file(name+".nobg")
        seq_graph = obg.SequenceGraph.from_file(name+".nobg.sequences")
        var_list = cls(indexed_interval, graph, seq_graph)
        var_list.create_lookup(name+"_variants.vcf")
        return var_list


class GenotypeMatrix:
    def __init__(self, matrix, positions):
        self._matrix = matrix
        self._positions = positions

    def find_path(self, genotype):
        variants = set(genotype.variants)
        start_idx = np.searchsorted(self._positions, genotype.start_pos)
        end_idx = np.searchsorted(self._positions, genotype.end_pos)
        print(start_idx, end_idx)
        relevant_matrix = self._matrix[start_idx:end_idx]
        relevant_ids = range(start_idx, end_idx)
        hits = np.array([variant_id in variants
                         for variant_id in relevant_ids])[:, None]
        return np.flatnonzero(np.all(relevant_matrix == hits, axis=0))

    @classmethod
    def from_vcf(cls, vcf_file_name):
        matrix = []
        positions = []
        for line in open(vcf_file_name):
            if line.startswith("#"):
                continue
            row = line.split("\t")
            positions.append(int(row[1]))
            vals = ["1" in val.split(":")[0] for val in vals[9:]]
            matrix.append(vals)
        return cls(np.array(matrix, dtype="bool"),
                   np.array(positions, dtype="int"))


def main(name, interval_file):
    var_list = VariantList.from_name(name)
    ineterval_collection = IntervalCollection.from_file(interval_file)
    types = [var_list.create_type(interval) for
             interval in interval_collection]
    matrix = GenoTypeMatrix.from_vcf("chr"+name+"_variants.vcf")
    hits_list = [matrix.find_path(genotype) for genotype in types]
    print(hits_list)
