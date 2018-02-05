import unittest
import math
import numpy as np
from graph_peak_caller.nongraphpeaks import NonGraphPeak, NonGraphPeakCollection


class TestNonGraphPeakCollection(unittest.TestCase):


    def setUp(self):
        # Create sample test bed file
        f = open("test_bed_file.bed", "w")
        f.writelines(["chr1\t1\t5\t.\t0\t+\t.\t.\t4.5\n"])
        f.writelines(["chr1\t10\t12\t.\t0\t+\t.\t.\t5.5\n"])
        f.close()

        # Create sample fasta file
        f = open("test_fasta_file.bed", "w")
        f.writelines([">chr1\n", "AAACCCTTTGGGAAACCCTTTGGG\n"])
        f.close()

    def test_from_bed_file(self):
        peaks = NonGraphPeakCollection.from_bed_file("test_bed_file.bed")
        peaks = peaks.peaks
        self.assertEqual(peaks[0].chromosome, "chr1")
        self.assertEqual(peaks[0].start, 1)
        self.assertEqual(peaks[0].score, 4.5)
        self.assertEqual(peaks[1].score, 5.5)
        self.assertEqual(peaks[1].end, 12)

    def test_filter_peaks_outside_region_other_chromosome(self):
        peaks = NonGraphPeakCollection.from_bed_file("test_bed_file.bed")
        peaks.filter_peaks_outside_region("chr2", 4, 7)
        self.assertEqual(len(peaks.peaks), 0)

    def test_filter_peaks_outside_region_same_chromosome(self):
        peaks = NonGraphPeakCollection.from_bed_file("test_bed_file.bed")
        peaks.filter_peaks_outside_region("chr1", 0, 11)
        self.assertEqual(len(peaks.peaks), 1)
        self.assertEqual(peaks.peaks[0].start, 1)

    def test_sort_on_sore(self):
        peaks = NonGraphPeakCollection.from_bed_file("test_bed_file.bed")
        peaks._sort_on_score_descending()
        self.assertTrue(np.isclose(peaks.peaks[0].score, 5.5))

    def test_set_sequences(self):
        peaks = NonGraphPeakCollection.from_bed_file("test_bed_file.bed")
        peaks.set_peak_sequences_using_fasta(fasta_file_location="test_fasta_file.bed")
        self.assertEqual(peaks.peaks[0].sequence, "AACC")
        self.assertEqual(peaks.peaks[1].sequence, "GG")

    def test_to_fasta(self):
        peaks = NonGraphPeakCollection.from_bed_file("test_bed_file.bed")
        peaks.set_peak_sequences_using_fasta("test_fasta_file.bed")
        peaks.save_to_sorted_fasta("test_fasta.fasta")

        f = open("test_fasta.fasta")
        lines = list(f.readlines())
        self.assertTrue(lines[0].startswith(">"))
        self.assertEqual(lines[3], "AACC\n")
        self.assertTrue(lines[2].startswith(">"))
        self.assertEqual(lines[1], "GG\n")

if __name__ == "__main__":
    unittest.main()
