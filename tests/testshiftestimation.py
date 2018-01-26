from collections import defaultdict
import logging
import numpy as np
import cProfile


class Treatment:
    def __init__(self, dicts):
        self._chroms = list(set(list(dicts["+"].keys()) + list(dicts["-"].keys())))
        self._pos_dict = {key: np.array(dicts["+"][key], dtype="int")
                          for key in self._chroms}
        self._neg_dict = {key: np.array(dicts["-"][key], dtype="int")
                          for key in self._chroms}
        for val in self._pos_dict.values():
            val.sort()
        for val in self._neg_dict.values():
            val.sort()
        self.total = sum(val.size for val in self._pos_dict.values())
        self.total += sum(val.size for val in self._neg_dict.values())

    def get_chr_names(self):
        return self._chroms
        return ["chr%s" % i for i in range(1, 10)]

    def get_locations_by_chr(self, chrom):
        return self._pos_dict[chrom], self._neg_dict[chrom]

        n_peaks = 10000
        a = np.array([1, 2, 3, 4, 5], dtype="int")
        pos_tags = np.empty(5*n_peaks, dtype="int")
        for i in range(n_peaks):
            pos_tags[(i*5):((i+1)*5)] = i*1000+a
        return pos_tags, pos_tags+50

    @classmethod
    def from_bed_file(cls, file_name):
        dicts = {"+": defaultdict(list),
                 "-": defaultdict(list)}
        with open(file_name) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                chr_name = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                strand = parts[5]
                pos = start if strand == "+" else end
                dicts[strand][chr_name].append(pos)

        return cls(dicts)


class Opt:
    def __init__(self):
        self.umfold = 50
        self.lmfold = 5
        self.bw = 300
        self.info = logging.info
        self.debug = logging.debug
        self.warn = logging.warning
        self.gsize = 300000000

t = Treatment.from_bed_file("vgdata/ctcf_reads_chr1.bed")
o = Opt()

# cProfile.run("p = PeakModel(o, t)")
# PeakModel(o, t)")
# print(p.d)
