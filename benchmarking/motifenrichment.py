from gendatafetcher.sequences import get_sequence_ucsc
from Bio import Seq, SeqIO, SeqRecord
from pybedtools import BedTool, Interval
import logging
logging.basicConfig(level=logging.INFO)
from nongraphpeaks import NonGraphPeak, NonGraphPeakCollection
import subprocess
import matplotlib.pyplot as plt


class MotifMatcher():

    def __init__(self, ranked_fasta_file_name, meme_motif_file_name):
        self.fasta_file = ranked_fasta_file_name
        self.meme_file = meme_motif_file_name

        self.peaks_matching_motif = set()
        self.sorted_peaks = []
        self.get_sorted_peak_ids()
        self.get_peaks_matching_motif()

    def get_peaks_matching_motif(self):
        #commmand = ["/home/ivargry/meme_4.11.4/src/fimo", self.meme_file, self.fasta_file]
        #commmand = ["bash", "-c", "'/home/ivargry/meme_4.11.4/src/fimo -oc fimo_tmp " + self.meme_file + " " + self.fasta_file +"'"]
        #print(' '.join(commmand))
        subprocess.check_output(["~/Sources/meme_4.12.0/src/fimo -oc fimo_tmp %s %s" % (self.meme_file, self.fasta_file)], shell=True)
        #ps = subprocess.check_output(commmand, shell=True)
        #output, error = ps.communicate()
        #print("Output")
        #print(output)
        #print("Error")
        #print(error)

        result_file = "fimo_tmp/fimo.txt"

        for line in open(result_file):
            l = line.split()
            sequence_id = l[2]
            self.peaks_matching_motif.add(sequence_id)

        #print(self.peaks_matching_motif)

    def get_sorted_peak_ids(self):
        for line in open(self.fasta_file):
            if line.startswith(">"):
                id = line.split()[0].replace(">", "")
                self.sorted_peaks.append(id)

    def compute_true_positives(self):
        true_positives = []
        n_matched = 0
        n_checked = 0
        for peak in self.sorted_peaks:
            if peak in self.peaks_matching_motif:
                n_matched += 1

            n_checked += 1
            true_positives.append(n_matched / n_checked)

        return true_positives


def plot_true_positives(peak_file_sets, meme_file_name):
    for name, fasta_file_name in peak_file_sets.items():
        matcher = MotifMatcher(fasta_file_name, meme_file_name)
        true_positives = matcher.compute_true_positives()
        n_matching = len(matcher.peaks_matching_motif)
        print("Fasta: %s" % fasta_file_name)
        n_tot = len(matcher.sorted_peaks)
        print("N tot: %d" % n_tot)
        print("True positives for %s: %d / %d = %.3f" % (name, n_matching, n_tot, n_matching/n_tot))
        plt.plot(true_positives, label=name + " (%.2f)" % (100 * n_matching/n_tot))
    plt.legend()
    plt.show()


if __name__ == "__main__":

    collection = NonGraphPeakCollection.from_bed_file("../tests/CTCF_peaks.narrowPeak")
    collection.filter_peaks_outside_region("chr6", 28510119, 33480577)
    collection.set_peak_sequences()
    collection.save_to_sorted_fasta("CTCF_filtered.fasta")
    #import sys
    #sys.exit()
    #matcher = MotifMatcher("../tests/real_data_sequences.fasta", "MA0139.1.meme")
    #true_positives = matcher.compute_true_positives()

    #print(true_positives)
    #sys.exit()
    plot_true_positives(
        {
            "graph peaks2": "../tests/real_data_sequences_q30.fasta",
            "graph peaks on macs path": "../tests/real_data_sequences_only_macs_path.fasta",
            "macs reads remapped": "../tests/real_data_sequences_macs_remapped.fasta",
            "macs": "CTCF_filtered.fasta",
            #"graph_peaks_bugfix": "../tests/real_data_sequences_after_bugfix"
            #"macs": "../tests/CTCF_peaks.narrowPeak"
            #"macs": "CTCFpFix.fasta"
            # "new_fasta": "../tests/sequences_new_control_sample.fasta",
            # "max_instead_of_average": "../tests/real_data_sequences_max_instead_of_average.fasta",
            #"spp_from_encode": "ENCFF155DHA_filtered.fasta"

        },
        "MA0139.1.meme"
    )
