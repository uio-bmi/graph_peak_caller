#from gendatafetcher.sequences import get_sequence_ucsc
from Bio import Seq, SeqIO, SeqRecord
from pybedtools import BedTool, Interval
import logging
from pyfaidx import Fasta
import json
from .peakcollection import PeakCollection
import numpy as np

logging.basicConfig(level=logging.INFO)


class NonGraphPeak():
    def __init__(self, chromosome, start, end, score=None, unique_id=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.score = score
        self.sequence = None
        self.unique_id = None

    def set_sequence_using_fasta_index(self, fasta_index):
        self.sequence = str(fasta_index[self.chromosome][self.start:self.end])

    def __str__(self):
        return "Peak(%s, %d, %d, score=%.3f)" % (self.chromosome,
                                                 self.start,
                                                 self.end,
                                                 self.score)

    def to_file_line(self):
        object = {
                "chromosome": self.chromosome,
                "start": int(self.start),
                "end": int(self.end),
                "score": float(self.score)
                  }
        try:
            d = json.dumps(object)
        except:
            for k, v in object.items():
                print(k, v, type(v))
            raise
        return d

    @classmethod
    def from_file_line(cls, line):
        object = json.loads(line)
        return cls(object["chromosome"], object["start"], object["end"],
                   score=object["score"])

    def __repr__(self):
        return self.__str__()


class NonGraphPeakCollection(object):
    def __init__(self, peaks):
        self.peaks = peaks

    @classmethod
    def from_bed_file(cls, file_name):
        peaks = []
        bedpeaks = BedTool(file_name)
        for peak in bedpeaks:
            score = float(peak[8])  # q value

            peaks.append(NonGraphPeak(peak.chrom,
                                      peak.start,
                                      peak.end,
                                      score))
        return cls(peaks)

    def get_graph_peak_collection(self, graph, linear_interval_through_graph, graph_region=None):
        return PeakCollection.create_from_nongraph_peak_collection(graph, self,
                                                                   linear_interval_through_graph,
                                                                   graph_region)


    def filter_peaks_outside_region(self, chromosome, start, end):
        new_peaks = []
        for peak in self.peaks:
            if peak.chromosome != chromosome:
                continue
            if peak.start < start or peak.end > end:
                continue
            new_peaks.append(peak)
        self.peaks = new_peaks

    def _sort_on_score_descending(self):
        self.peaks.sort(key=lambda peak: peak.score, reverse=True)

    def to_fasta(self, file_name):
        """
        lines = (SeqRecord.SeqRecord(Seq.Seq(peak.sequence),
                                     id=str(i), description=peak.to_file_line())
                 for i, peak in enumerate(self.peaks))

        SeqIO.write(lines, file_name, "fasta")
        """
        f = open(file_name, "w")
        for i, peak in enumerate(self.peaks):
            if peak.unique_id is None:
                id = i
            else:
                id = peak.unique_id
            f.writelines([">%d %s\n" % (id, peak.to_file_line())])
            f.writelines(["%s\n" % peak.sequence])
        f.close()

    @classmethod
    def from_fasta(cls, file_name):
        f = open(file_name)
        peaks = []
        while True:
            header = f.readline()
            sequence = f.readline()
            if not sequence:
                break

            header = header.split(maxsplit=1)
            id = header[0].replace(">", "")
            interval_json = header[1]
            peak = NonGraphPeak.from_file_line(interval_json)
            peak.unique_id = id
            assert peak.unique_id is not None
            peak.sequence = sequence
            peaks.append(peak)

        avg_peak_size = np.mean([p.end - p.start for p in peaks])
        logging.info("Avg peak size: %.2f" % avg_peak_size)

        return cls(peaks)

    def set_peak_sequences_using_fasta(self, fasta_file_location="grch38.fasta"):
        logging.info("Setting peak sequences using fasta index")
        genome = Fasta(fasta_file_location)
        i = 0
        for peak in self.peaks:
            if i % 10000 == 0:
                logging.info("%d/%d peaks processed" % (i, len(self.peaks)))
            i += 1

            peak.set_sequence_using_fasta_index(genome)

    def save_to_sorted_fasta(self, file_name):
        self._sort_on_score_descending()
        self.to_fasta(file_name)
