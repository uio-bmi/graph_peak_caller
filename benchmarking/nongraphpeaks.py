from gendatafetcher.sequences import get_sequence_ucsc
from Bio import Seq, SeqIO, SeqRecord
from pybedtools import BedTool, Interval
import logging
logging.basicConfig(level=logging.INFO)


class NonGraphPeak():
    def __init__(self, chromosome, start, end, score=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.score = score
        self.sequence = None

    def set_sequence(self):
        self.sequence = get_sequence_ucsc(self.chromosome, self.start, self.end)

    def __str__(self):
        return "Peak(%s, %d, %d, score=%.3f)" % (self.chromosome,
                                                 self.start,
                                                 self.end,
                                                 self.score)

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
        lines = (SeqRecord.SeqRecord(Seq.Seq(peak.sequence),
                 id=str(i), description=str(peak))
                 for i, peak in enumerate(self.peaks))

        SeqIO.write(lines, file_name, "fasta")

    def set_peak_sequences(self):
        i = 0
        for peak in self.peaks:
            logging.info("Set sequence for peak %d/%d" % (i, len(self.peaks)))
            i += 1
            peak.set_sequence()

    def save_to_sorted_fasta(self, file_name):
        self._sort_on_score_descending()
        self.to_fasta(file_name)
