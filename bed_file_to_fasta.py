from pybedtools import BedTool, Interval
from gendatafetcher.sequences import get_sequence_ucsc
from Bio import Seq, SeqIO, SeqRecord

MHC_REGION = Interval("chr6", 28510119, 33480577)


def bed_to_fasta(filename, outfile_name=None):
    if outfile_name is None:
        outfile_name = filename + ".fa"
    peaks = BedTool(filename)
    peaks = peaks.filter(lambda x: x.chrom == MHC_REGION.chrom and
                         x.start >= MHC_REGION.start and
                         x.end <= MHC_REGION.end)
    sequences = (get_sequence_ucsc(peak.chrom, peak.start, peak.end, False)
                 for peak in peaks)
    sequences = (SeqRecord.SeqRecord(Seq.Seq(seq), id=str(i), description="peak")
                 for i, seq in enumerate(sequences))
    SeqIO.write(sequences, outfile_name, "fasta")


if __name__ == "__main__":
    bed_to_fasta("/home/knut/Data/ENCODE_peaks/MACS/ENCSR000AHD_MCF-7_CTCF_peaks.narrowPeak",
                 "test_fasta.fa")
