import numpy as np
from pyfaidx import Fasta
import offsetbasedgraph as obg
import logging
from collections import namedtuple, deque
from itertools import chain


Variant = namedtuple("Variant", ["offset", "ref", "alt"])
FullVariant = namedtuple("FullVariant", ["offset", "ref", "alt", "precence"])


class VariantPrecence:
    def __init__(self, precence):
        self._precence = precence

    def get_samples(self, variant, f=None):
        s = self._precence if f is None else self._precence[f, :]
        match = np.any(s == variant, axis=-1)
        return np.flatnonzero(match)

    @classmethod
    def from_line(cls, line):
        parts = line.replace(".", "0").replace("/", "|").split("\t")
        haplo_type = [part.split(":", 2)[0].split("|") for part in parts]
        try: 
            precence = np.array([[int(a), int(b)] for a, b in haplo_type])
        except:
            print(line)
            print(haplo_type)
            raise
        return cls(precence)

    @classmethod
    def join(cls, precences, counts):
        offsets = [0] + list(np.cumsum(counts))[:-1]
        combined = np.zeros_like(precences[0]._precence)
        for precence, offset in zip(precences, offsets):
            combined[precence._precence > 0] += precence._precence[precence._precence > 0] + offset
        return cls(combined)


class VariantList:
    def __init__(self, start):
        self._start = int(start)
        self._variants = []

    def append(self, variant):
        self._variants.append(variant)

    def finalize(self):
        final_list = []
        cur_end = 0
        var_buffer = []
        for variant in self._variants:
            transformed_var = FullVariant(variant.offset-self._start,
                                      variant.ref, variant.alt,
                                      variant.precence)
            if transformed_var.offset >= cur_end:
                if var_buffer:
                    final_list.append(self._join_variants(var_buffer))
                var_buffer = []
            var_buffer.append(transformed_var)
            cur_end = max(transformed_var.offset + len(transformed_var.ref), cur_end)
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

    def get_variants_from_intervals(self, intervals):
        """ intervals reveresely sorted on start """
        intervals = iter(intervals)
        current_intervals = deque([])
        cur_interval = next(intervals)
        is_finished = False
        for line in self.f:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 9)
            try:
                pos = int(parts[1])-1
            except:
                continue
            while pos >= cur_interval[0] and not is_finished:
                current_intervals.append([cur_interval[-1], False,
                                          VariantList(cur_interval[0])])
                try:
                    cur_interval = next(intervals)
                except StopIteration:
                    is_finished = True
            if not current_intervals:
                continue

            # end = pos + len(parts[3])
            for entry in current_intervals:
                entry[1] = pos >= entry[0]

            while current_intervals and current_intervals[0][1]:
                yield current_intervals.popleft()[-1].finalize()
            if not current_intervals:
                if is_finished:
                    break
                continue
            ref = parts[3]
            alt = parts[4].lower().split(",")
            precence = VariantPrecence.from_line(parts[-1])
            var = FullVariant(int(pos), ref.lower(), alt, precence)
            for _, _, entry in current_intervals:
                entry.append(var)

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
    if ref_seq == alt_seq:
        return ["REF"]
    tentative_valid = [([], 0, 0)]
    for j, variant in enumerate(variants):
        next_tentative = []
        logging.debug(variant)
        for tentative in tentative_valid:
            cur_vars, alt_offset, prev_offset = tentative
            logging.debug(tentative)
            if not alt_seq[prev_offset+alt_offset:variant.offset+alt_offset] == ref_seq[prev_offset:variant.offset]:
                logging.debug("--------")
                continue
            for code, seq in enumerate([variant.ref] + variant.alt):
                logging.debug(alt_seq[variant.offset+alt_offset:variant.offset+alt_offset+len(seq)], seq.lower())
                if alt_seq[variant.offset+alt_offset:variant.offset+alt_offset+len(seq)] == seq.lower():
                    logging.debug("!")
                    next_tentative.append((cur_vars+[code], alt_offset+len(seq)-len(variant.ref), variant.offset+len(variant.ref)))
        tentative_valid = next_tentative
        assert all(len(t[0]) == j+1 for t in tentative_valid), (j, tentative_valid)
    real = []
    for tentative in tentative_valid:
        cur_vars, alt_offset, prev_offset = tentative
        if not len(alt_seq)-alt_offset == len(ref_seq):
            continue
        if not alt_seq[prev_offset+alt_offset:len(ref_seq)+alt_offset] == ref_seq[prev_offset:len(ref_seq)]:
            continue
        real.append(tentative)

    if not len(real) == 1:
        logging.debug(alt_seq)
        logging.debug(ref_seq)
        logging.debug(variants)
        logging.debug("---->")
        logging.debug(real, "#", tentative_valid)
        return ["No Match"]
    codes = real[0][0]
    f = None
    for code, variant in zip(codes, variants):
        f = variant.precence.get_samples(code, f)
    return f


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

    def run_peaks(self, peaks):
        alts, refs, intervals = zip(*(self.get_sequence_pair(peak) for peak in peaks))
        haplotypes = [traverse_variants(alt.lower(), ref.lower(), variants) for alt, ref, variants in
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

    def get_sequence_pair(self, peak):
        peak.graph = self.graph
        seq = self.seq_graph.get_interval_sequence(peak)
        assert len(seq) == peak.length(), (len(seq), peak.length())
        seq, interval = self.extend_peak(seq, peak.start_position, peak.end_position)
        ref_seq = self.fasta[int(interval[0]):int(interval[1])].seq
        if all(rp in self.indexed_interval.nodes_in_interval() for rp in peak.region_paths):
            if len(seq) == len(ref_seq):
                pass  # assert seq == ref_seq.lower(), (seq, ref_seq.lower())
        return seq, ref_seq, interval

    def extend_peak(self, seq, start_position, end_position):
        seq, start_offset = self.extend_peak_start(
            seq, start_position.offset, start_position.region_path_id)
        seq, end_offset = self.extend_peak_end(
            seq, end_position.offset, end_position.region_path_id)
        return seq, (start_offset, end_offset)

    def extend_peak_start(self, seq, start_offset, start_node):
        if start_node in self.indexed_interval.nodes_in_interval():
            return seq, self.indexed_interval.get_offset_at_node(start_node) + start_offset
        prev_node = max(-node for node in self.graph.reverse_adj_list[-start_node]
                        if -node in self.indexed_interval.nodes_in_interval())
        p_len = self.graph.node_size(prev_node)
        seq = self.seq_graph.get_interval_sequence(obg.Interval(p_len-1, start_offset, [prev_node, start_node]))+seq
        start = self.indexed_interval.get_offset_at_node(prev_node)+p_len-1
        return seq, start

    def extend_peak_end(self, seq, end_offset, end_node):
        if end_node in self.indexed_interval.nodes_in_interval():
            return seq, self.indexed_interval.get_offset_at_node(end_node) + end_offset
        next_node = min(node for node in self.graph.adj_list[end_node]
                        if node in self.indexed_interval.nodes_in_interval())
        p_len = self.graph.node_size(next_node)
        seq += self.seq_graph.get_interval_sequence(obg.Interval(end_offset, 1, [end_node, next_node]))
        end = self.indexed_interval.get_offset_at_node(next_node)+1
        return seq, end

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
    h = "./.:.:.	1|1:40:8	1|1:40:11	0|0:40:9	./.:.:.	./.:.:.	./.:.:.	0|0:40:8	0|0:40:9	./.:.:.	0|0:40:12	0|0:40:16	0|0:40:40	0|0:40:35	0|0:40:21	0|0:40:13	0|0:40:44	0|0:40:15	0|0:40:86	0|0:40:69	./.:.:.	0|0:40:14	0|0:40:59	0|0:40:51"
    precence = VariantPrecence.from_line(h)
    print(precence._precence)
    print(precence.get_samples(1, [0, 1, 3]))

    #alt = "aaaaaataagacgt"
    #ref = "ataaataagacgt"
    #variants = [Variant(offset=1, ref='T', alt=['A']), Variant(offset=8, ref='GACGTACCCTCA', alt=['G', 'GAGGTACCCTCA'])]

    # alt = "cccttctttttttg"
    # ref = "cttttttttttttg"
    # variants = [Variant(offset=0, ref='CTT', alt=['TTT', 'C', 'CTTTT', 'CTTC']), Variant(offset=5, ref='T', alt=['C']), Variant(offset=6, ref='T', alt=['A']), Variant(offset=7, ref='T', alt=['C', 'TC'])]

    # alt = "ggaaataaaaaa"
    # ref = "ggaaataaaaa"
    # variants = [Variant(offset=5, ref='T', alt=['C', 'TA']), Variant(offset=7, ref='A', alt=['C']), Variant(offset=8, ref='A', alt=['T', 'G']), Variant(offset=11, ref='T', alt=['A'])]

    # alt = "agtttcactggg"
    # ref = "agtttcagtagg"
    # variants = [Variant(offset=7, ref='GTA', alt=['CTA', 'ATA', 'GT', 'GTG'])]
    # [Variant(offset=7, ref='G', alt=['C', 'A']), Variant(offset=8, ref='TA', alt=['T']), Variant(offset=9, ref='A', alt=['G'])]

#     alt = "accttatagaaa"
#     ref = "accttataagaaa"
#     variants = [Variant(offset=6, ref='TA', alt=['T'])]

    # alt = "agtttcactggg"
    # ref = "agtttcagtagg"
    # variants = [Variant(offset=7, ref='G', alt=['C', 'A']), Variant(offset=9, ref='A', alt=['', 'G'])]
    # print(traverse_variants(alt, ref, variants))
