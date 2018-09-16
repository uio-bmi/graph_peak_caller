import numpy as np
from pyfaidx import Fasta
import offsetbasedgraph as obg
import logging
from collections import namedtuple, deque
from itertools import chain, dropwhile, takewhile


Variant = namedtuple("Variant", ["offset", "ref", "alt"])
FullVariant = namedtuple("FullVariant", ["offset", "ref", "alt", "precence"])


def find_valid_haplotypes(variants, combination):
    combination = np.asanyarray(combination)
    if not np.any(combination):
        return np.arange(1136)
    nonzero = np.flatnonzero(combination)
    informative_haplotypes = []
    informative_variants = [variant for i, variant in enumerate(variants)
                            if i in nonzero]
    for i, variant in enumerate(informative_variants):
        informative_haplotypes.extend(
            variant.precence.get_samples(combination[nonzero[i]]))
    f = np.unique(informative_haplotypes)
    if not f.size:
        return np.array([])
    for code, variant in zip(combination, variants):
        f = variant.precence.get_samples(code, f)
    return f


class VariantPrecence:
    strict = False
    accept_ref = False

    def __init__(self, precence):
        self._precence = precence

    def get_samples(self, variant, f=None):
        s = self._precence if f is None else self._precence[f, :]
        if self.accept_ref and variant == 0:
            return f if f is not None else np.arange(self._precence.shape[0])
        match = np.any((s == variant) | (s == -1),
                       axis=-1)
        res = np.flatnonzero(match)
        if f is not None:
            return f[res]
        return res

    @classmethod
    def from_line(cls, line):
        if cls.strict:
            parts = line.replace(".", "0").replace("/", "|").split("\t")
        else:
            parts = line.replace(".", "-1").replace("/", "|").split("\t")
        haplo_type = [part.split(":", 2)[0].split("|") for part in parts]
        precence = np.array([[int(a), int(b)] for a, b in haplo_type])
        return cls(precence)

    @classmethod
    def join(cls, precences, counts):
        offsets = [0] + list(np.cumsum(counts))[:-1]
        combined = np.zeros_like(precences[0]._precence)
        for precence in precences:
            combined[precence._precence < 0] = -1
        for precence, offset in zip(precences, offsets):
            combined[precence._precence > 0] = precence._precence[precence._precence > 0] + offset
        return cls(combined)

    def __repr__(self):
        return "P(%s)" % np.count_nonzero(self._precence > 0)


class VariantList:
    def __init__(self, start, end):
        self._start = int(start)
        self._end = int(end)
        self._variants = []

    def append(self, variant):
        self._variants.append(variant)

    def finalize(self):
        return [FullVariant(variant.offset-self._start,
                            variant.ref, variant.alt,
                            variant.precence) for variant in self._variants]

        final_list = []
        cur_end = 0
        var_buffer = []

        for variant in self._variants:
            transformed_var = FullVariant(
                variant.offset-self._start,
                variant.ref, variant.alt,
                variant.precence)
            assert transformed_var.offset <= self._end-self._start, (variant, self._start, transformed_var.offset, self._end-self._start)
            if transformed_var.offset >= cur_end:
                if var_buffer:
                    final_list.append(self._join_variants(var_buffer))
                var_buffer = []
            var_buffer.append(transformed_var)
            cur_end = max(transformed_var.offset + len(transformed_var.ref),
                          cur_end)
        if var_buffer:
            final_list.append(self._join_variants(var_buffer))
        return final_list

    def _join_variants(self, variants):
        if len(variants) == 1:
            return variants[0]
        start = min(variant.offset for variant in variants)
        end = max(variant.offset + len(variant.ref) for variant in variants)
        ref = []*(end-start)
        for variant in variants:
            ref[variant.offset-start:variant.offset-start+len(variant.ref)] = variant.ref
        ref = "".join(ref)
        alts = [ref[:variant.offset-start] + v_alt + ref[variant.offset+len(variant.ref)-start:]
                for variant in variants for v_alt in variant.alt]
        counts = [len(variant.alt) for variant in variants]
        precences = [variant.precence for variant in variants]
        combined_precence = VariantPrecence.join(precences, counts)
        return FullVariant(start, ref, alts, combined_precence)
def prune_seqs(ref, alts):
    offset = 0
    for cs in zip(*([ref]+alts)):
        if all(c == cs[0] for c in cs):
            offset += 1
        else:
            break
    return ref[offset:], [alt[offset:] for alt in alts], offset


class VCF:
    def __init__(self, vcf_file_name):
        self.file_name = vcf_file_name
        self.f = open(self.file_name)
        self.cur_lines = []

    def _join_variants(self, variants):
        if len(variants) == 1:
            return variants[0]
        start = min(variant.offset for variant in variants)
        end = max(variant.offset + len(variant.ref) for variant in variants)
        ref = []*(end-start)
        for variant in variants:
            ref[variant.offset-start:variant.offset-start+len(variant.ref)] = variant.ref
        ref = "".join(ref)
        alts = [ref[:variant.offset-start] + v_alt + ref[variant.offset+len(variant.ref)-start:]
                for variant in variants for v_alt in variant.alt]
        return Variant(start, ref, alts)


    def get_haplotype_sequences(self, haplotype, interval, ref_seq, lines=None):
        interval = (int(interval[0]), int(interval[1]))
        if lines is None:
            lines = self.f
        line_parts = (line.split("\t", 9) for line in lines
                      if not line.startswith("#"))
        interval_line_parts = dropwhile(
            lambda parts: int(parts[1])-1 < interval[0],
            line_parts)
        interval_line_parts = takewhile(
            lambda parts: int(parts[1])-1 <= interval[1],
            interval_line_parts)

        seqs = [[], []]
        cur_positions = [0, 0]
        for parts in interval_line_parts:
            variant = VariantPrecence.from_line(parts[-1])._precence[haplotype]
            variant_start = int(parts[1])-1-int(interval[0])
            variant_end = variant_start + len(parts[3])
            alts = parts[4].split(",")
            for i, (seq, allele) in enumerate(zip(seqs, variant)):
                if allele > 0:
                    seq.append(ref_seq[cur_positions[i]:variant_start])
                    seq.append(alts[allele-1])
                    cur_positions[i] = variant_end
        for pos, seq in zip(cur_positions, seqs):
            seq.append(ref_seq[pos:interval[1]-interval[0]])
        return ["".join(seq) for seq in seqs]

    def get_variants_from_intervals(self, intervals):
        """ intervals reveresely sorted on start """
        PRUNE = True
        StackElement = namedtuple("StackElement", ["end", "variant_list"])
        intervals = iter(intervals)
        current_intervals = deque([])
        cur_interval = next(intervals)
        is_finished = False
        line_parts = (line.split("\t", 9) for line in self.f
                      if not line.startswith("#"))
        for parts in line_parts:
            pos = int(parts[1])-1
            while pos >= cur_interval[0] and not is_finished:
                current_intervals.append(
                    StackElement(cur_interval[-1],
                                 VariantList(*cur_interval)))
                try:
                    cur_interval = next(intervals)
                except StopIteration:
                    is_finished = True
            if not current_intervals:
                continue

            while current_intervals and current_intervals[0].end < pos:
                tmp = current_intervals.popleft()
                yield tmp.variant_list.finalize()
            if not current_intervals:
                if is_finished:
                    break
                continue
            ref = parts[3].lower()
            alt = parts[4].lower().split(",")
            if PRUNE:
                ref, alt, offset = prune_seqs(ref, alt)
                pos = int(pos) + offset
            precence = VariantPrecence.from_line(parts[-1])
            var = FullVariant(int(pos), ref.lower(), alt, precence)
            variant_end = pos+len(ref)
            for stack_element in current_intervals:
                if stack_element.end > pos:
                    stack_element.variant_list.append(var)

    def get_variants_from(self, start, end):
        variants = []
        cur_end = 0
        intersecting = []
        for line in chain(self.cur_lines, self.f):
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            pos = int(parts[1])-1
            if pos < start:
                continue
            if pos >= cur_end:
                if intersecting:
                    variants.append(self._join_variants(intersecting))
                intersecting = []
            if pos >= end:
                self.cur_lines = [line]
                break

            ref = parts[3]
            intersecting.append(Variant(int(pos-start), ref, parts[4].split(",")))
            cur_end = max(pos+len(ref), cur_end)
            continue
        if intersecting:
            variants.append(self._join_variants(intersecting))

        self.cur_lines = []
        variants.sort(key=lambda x: x.offset)
        return variants


def traverse_variants(alt_seq, ref_seq, variants):
    if not len(variants):
        return None
    Combination = namedtuple("Combination", ["code", "alt_offset", "prev_offset"])
    combinations = [Combination([], 0, 0)]
    for j, variant in enumerate(variants):
        new_combinations = []
        for comb in combinations:
            inter_alt = alt_seq[comb.prev_offset+comb.alt_offset:variant.offset+comb.alt_offset]
            inter_ref = ref_seq[comb.prev_offset:variant.offset]
            if not inter_alt == inter_ref:
                continue
            new_combinations.append(Combination(
                comb.code+[0],
                comb.alt_offset,
                comb.prev_offset))
            if comb.prev_offset > variant.offset:
                continue
            for code, seq in enumerate(variant.alt):
                p_ref, p_alt, p_offset = prune_seqs(variant.ref, [seq])
                seq = p_alt[0]
                start_on_alt = variant.offset+comb.alt_offset+p_offset
                var_alt = alt_seq[start_on_alt:start_on_alt+len(seq)]
                if var_alt == seq.lower():
                    new_combinations.append(Combination(
                        comb.code+[code+1],
                        comb.alt_offset+len(seq)-len(p_ref),
                        variant.offset+len(variant.ref)))

        combinations = new_combinations
        # assert all(len(t[0]) == j+1 for t in combinations), (j, combinations)
    real = []
    for comb in combinations:
        # if not len(alt_seq)-comb.alt_offset == len(ref_seq):
        #     continue
        alt_stub = alt_seq[comb.prev_offset+comb.alt_offset:len(ref_seq)+comb.alt_offset]
        ref_stub = ref_seq[comb.prev_offset:len(ref_seq)]
        l = min(len(alt_stub), len(ref_stub))
        if not alt_stub[:l] == ref_stub[:l]:
            continue
        real.append(comb)

    if not len(real):
        logging.info(alt_seq)
        logging.info(ref_seq)
        logging.info(variants)
        logging.info("---->")
        logging.info("%s / %s", real, combinations)
        print("")
        return ["No Match"]

    haplotypes = []
    for valid in real:
        codes = valid[0]
        haplotypes.extend(find_valid_haplotypes(variants, codes))
        # f = None
        # for code, variant in zip(codes, variants):
        #     f = variant.precence.get_samples(code, f)
        # haplotypes.extend(f)

    return np.sort(np.unique(haplotypes))
    # codes = real[0][0]
    # f = None
    # for code, variant in zip(codes, variants):
    #     f = variant.precence.get_samples(code, f)
    # return f


def find_haplotype(seq, refseq, vcf, start, end):
    refseq = refseq.lower()
    if seq == refseq:
        return ["REF"]
    # logging.info("Finding haplotype for %s->%s, (%s, %s)", refseq, seq, start, end)
    return traverse_variants(seq, refseq, vcf.get_variants_from(start, end))
    cur_offset = 0
    cur_ref_offset = 0
    haplotypes = []

    possible_vars = []
    for i, variant in enumerate(vcf.get_variants_from(start, end)):
        possible_vars.append(variant)
        new_offset = variant.offset
        l = new_offset - cur_ref_offset
        if seq[cur_offset:cur_offset+l] != refseq[cur_ref_offset:new_offset]:
            logging.error("%s != %s", seq[cur_offset:cur_offset+l], refseq[cur_ref_offset:new_offset])
            return ["ERROR"]
        cur_offset += l
        cur_ref_offset = new_offset
        if seq[cur_offset:cur_offset+len(variant.alt)] == variant.alt.lower():
            if (not seq[cur_offset:cur_offset+len(variant.ref)] == variant.ref.lower()) or len(variant.alt) > len(variant.ref):
                haplotypes.append(1)
                cur_offset += len(variant.alt)
                cur_ref_offset += len(variant.ref)
                continue
        # assert seq[cur_offset:cur_offset+len(variant.ref)] == variant.ref
        haplotypes.append(0)
    logging.info(possible_vars)
    logging.info([var for var, h in zip(possible_vars, haplotypes) if h])
    return haplotypes


class Main:

    def __init__(self, indexed_interval, graph, seq_graph, fasta, vcf):
        self.indexed_interval = indexed_interval
        self.graph = graph
        self.seq_graph = seq_graph
        self.fasta = fasta
        self.vcf = vcf
        self.counter = 0
        self.prev_end = -1

    def get_haplo_sequence_around_intervals(self, intervals, haplotype):
        extended = (self.extend_peak_to_reference(interval) for interval in intervals)
        ref_intervals = [
            (self.indexed_interval.get_offset_at_node(peak.region_paths[0]) + peak.start_position.offset,
             self.indexed_interval.get_offset_at_node(peak.region_paths[-1]) + peak.end_position.offset)
            for peak in extended]
        starts, ends = zip(*ref_intervals)
        interval = (min(starts), max(ends))
        ref_seq = self.fasta[int(interval[0]):int(interval[1])].seq
        var_seqs = self.vcf.get_haplotype_sequences(haplotype, interval, ref_seq)
        return var_seqs

    def run_peak_set(self, peaks_dict):
        all_name_tuples = [(key, self.get_sequence_pair(read))
                           for key, reads in peaks_dict.items() for read in reads]
        names, tuples = zip(*sorted(all_name_tuples, key=lambda x: x[1][2][0]))
        alts, refs, intervals = zip(*tuples)
        params = zip(alts, refs, self.vcf.get_variants_from_intervals(intervals))
        haplotypes = (traverse_variants(*param) for param in params)
        haplotypes = (h if h is not None else list(range(1136)) for h in haplotypes)
        result_dict = {name: [] for name in peaks_dict}
        for name, result in zip(names, haplotypes):
            
            result_dict[name].append(result)
        return result_dict

    def run_peaks(self, peaks):
        alts, refs, intervals = zip(*(self.get_sequence_pair(peak) for peak in peaks))
        haplotypes = [traverse_variants(alt.lower(), ref.lower(), variants)
                      for alt, ref, variants in
                      zip(alts, refs, self.vcf.get_variants_from_intervals(intervals))]
        return haplotypes

    def run_peak(self, peak):
        seq, ref_seq, interval = self.get_sequence_pair(peak)
        if interval[0] < self.prev_end:
            return []
        self.prev_end = interval[-1]
        haplo = find_haplotype(seq, ref_seq, self.vcf, interval[0], interval[1])
        if haplo:
            self.counter += 1
        return haplo

    def extend_peak_to_reference(self, peak):
        if peak.region_paths[0] < 0:
            peak.graph = self.graph
            peak = peak.get_reverse()
        start_node = peak.region_paths[0]
        end_node = peak.region_paths[-1]
        start_offset = peak.start_position.offset
        end_offset = peak.end_position.offset
        rps = peak.region_paths
        if start_node not in self.indexed_interval.nodes_in_interval():
            start_node = min(
                -node for node in self.graph.reverse_adj_list[-start_node]
                if -node in self.indexed_interval.nodes_in_interval())
            start_offset = self.graph.node_size(start_node)-1
            rps = [start_node] + rps
        if end_node not in self.indexed_interval.nodes_in_interval():
            end_node = max(node for node in self.graph.adj_list[end_node]
                           if node in self.indexed_interval.nodes_in_interval())
            end_offset = 1
            rps = rps + [end_node]
        return obg.Interval(start_offset, end_offset, rps)

    def get_sequence_pair(self, peak):
        peak.graph = self.graph
        peak = self.extend_peak_to_reference(peak)
        seq = self.seq_graph.get_interval_sequence(peak)
        assert len(seq) == peak.length(), (len(seq), peak.length())
        interval = (self.indexed_interval.get_offset_at_node(peak.region_paths[0]) + peak.start_position.offset,
                    self.indexed_interval.get_offset_at_node(peak.region_paths[-1]) + peak.end_position.offset)
        ref_seq = self.fasta[int(interval[0]):int(interval[1])].seq
        assert len(ref_seq) == interval[1]-interval[0]
        return seq.lower(), ref_seq.lower(), interval

    @classmethod
    def from_name(cls, name, fasta_file_name, chrom):
        indexed_interval = obg.NumpyIndexedInterval.from_file(
            name+"_linear_pathv2.interval")
        graph = obg.Graph.from_file(name+".nobg")
        seq_graph = obg.SequenceGraph.from_file(name+".nobg.sequences")
        fasta = Fasta(fasta_file_name)[chrom]
        vcf = VCF(name + "_variants.vcf")
        return cls(indexed_interval, graph, seq_graph, fasta, vcf)
    # var_list.create_lookup(name+"_variants.vcf")
    # return var_list


class DummyVCF:

    @staticmethod
    def get_variants_from(start):
        return [Variant(2, "A", "T"), Variant(4, "G", ""), Variant(5, "G", "C")]


def test():
    refseq = "AAAGGGTTT"
    seq = "AATGGCTTT"
    assert find_haplotype(seq, refseq, DummyVCF, 0) == [1, 0, 1]
    print("SUCCESS")

def test2():
    refseq = "AAAGGGTTT"
    seq = "AAAGGCTTT"
    assert find_haplotype(seq, refseq, DummyVCF, 0) == [0, 0, 1]
    print("SUCCESS")

def test3():
    refseq = "AAAGGGTTT"
    seq = "AAAGCTTT"
    assert find_haplotype(seq, refseq, DummyVCF, 0) == [0, 1, 1]
    print("SUCCESS")


def test4():
    vcf = "/home/knut/Documents/phd/data/graphs/1_variants.vcf"
    vcf = VCF(vcf)
    for i, var in enumerate(vcf.get_variants_from(70)):
        print(var)
        if i > 10:
            break


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    line = ["_", "10", "_",  "aaa", "ccc"] + ["_"]*4
    h = "./.:.:.	1|1:40:8	1|1:40:11	0|0:40:9	./.:.:.	./.:.:.	./.:.:.	0|0:40:8	0|0:40:9	./.:.:.	0|0:40:12	0|0:40:16	0|0:40:40	0|0:40:35	0|0:40:21	0|0:40:13	0|0:40:44	0|0:40:15	0|0:40:86	0|0:40:69	./.:.:.	0|0:40:14	0|0:40:59	0|0:40:51"
    line = "\t".join(line+ [h])
    print(line)
    ref = "ttaaaggcc"
    interval = (7, 7+len(ref))
    print(VCF("tmp.txt").get_haplotype_sequences(1, interval, ref, [line]))
#     
# 
#     alt = "ccttctttttttg"
#     ref = "ccttctttttttg"
#     variants = [FullVariant(offset=0, ref='ctt', alt=['ttt', 'c', 'ctttt'], precence=None),
#                 FullVariant(offset=2, ref='t', alt=['tc'], precence=None),
#                 FullVariant(offset=5, ref='t', alt=['c'], precence=None),
#                 FullVariant(offset=6, ref='t', alt=['a'], precence=None),
#                 FullVariant(offset=7, ref='t', alt=['c', 'tc'], precence=None)]
#     print(traverse_variants(alt, ref, variants))
