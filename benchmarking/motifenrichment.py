import logging
logging.basicConfig(level=logging.INFO)
#from nongraphpeaks import NonGraphPeak, NonGraphPeakCollection
import subprocess
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


class MotifMatcher():

    def __init__(self, ranked_fasta_file_name, meme_motif_file_name, run_fimo=True):
        self.fasta_file = ranked_fasta_file_name
        self.meme_file = meme_motif_file_name
        self.run_fimo = run_fimo
        self.peaks_matching_motif = set()
        self.sorted_peaks = []

        self.get_sorted_peak_ids()

        self.get_peaks_matching_motif()

    def get_peaks_matching_motif(self):

        out_dir = "fimo_" + self.fasta_file.replace(".fasta", "")

        if self.run_fimo:
            subprocess.check_output(["fimo -oc %s %s %s" % (out_dir, self.meme_file, self.fasta_file)], shell=True)

        result_file = out_dir + "/fimo.txt"

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
            else:
                print("Not match: %s" % peak)

            n_checked += 1
            true_positives.append(n_matched / n_checked)

        return true_positives


def plot_true_positives(peak_file_sets, meme_file_name, save_to_file=None, run_fimo=True):
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
    if save_to_file is not None:
        plt.savefig(save_to_file)
        print("Saved figure to %s" % save_to_file)
    else:
        plt.show()



if __name__ == "__main__":
    """
    from nongraphpeaks import NonGraphPeak, NonGraphPeakCollection
    collection = NonGraphPeakCollection.from_bed_file("../tests/macs_without_control_peaks.narrowPeak")
    collection.filter_peaks_outside_region("chr6", 28510119, 33480577)
    #collection.filter_peaks_outside_region("chr19", 54025634, 55084318)
    collection.set_peak_sequences()
    collection.save_to_sorted_fasta("../tests/mhc/macs_without_control.fasta")
    exit()
    """

    #import sys
    #sys.exit()
    #matcher = MotifMatcher("../tests/real_data_sequences.fasta", "MA0139.1.meme")
    #true_positives = matcher.compute_true_positives()

    #print(true_positives)
    #sys.exit()
    """
    plot_true_positives({
        #"macs with control": "../tests/macs_with_control_sequences_chr1.fasta",
        #"graph with control": "../tests/chr1_ctcf_with_control_sequences.fasta",
        #"graph": "../tests/chr1_ctcf_sequences.fasta",
        "graph mapped to whole genome": "../tests/sequences_mapped_whole_genome.fasta",
        #"graph using trimmed sequences": "../tests/using_trimmed_sequences.fasta",
        "macs": "../tests/macs_chr1.fasta"
    },
    "MA0139.1.meme")
    exit()
    """

    # SRF
    """
    plot_true_positives(
        {
            "macs": "macs_srf_without_control.fasta",
            "graph": "../tests/srf_sequences.fasta"
        },
        "MA0083.3.meme"
    )
    exit()
    """
    """
    plot_true_positives({
        #"graph without control": "../tests/mhc/ctcf_without_control_sequences.fasta",
        "graph with control": "../tests/ctcf_q50_with_control_sequences.fasta",
        "graph without control, stricter lambda": "../tests/mhc/ctcf_without_control_hack_sequences.fasta",
        #"old graph without control": "../tests/ctcf_q50_without_control_sequences.fasta",
        "macs without control": "macs_no_control_peaks.fasta",
        "graph using macs reads old": "../tests/linear_reads_moved_to_graph_sequences.fasta",
        "graph using macs reads": "../tests/mhc_using_macs_reads_sequences.fasta",
        "graph using macs reads remapepd": "../tests/mhc/macs_reads_remapped_sequences.fasta"
    },
    "MA0139.1.meme")
    """

    plot_true_positives(
        {
            "macs": "../tests/lrc_kir/macs_without_control.fasta",
            "graph without control": "../tests/lrc_kir/test_sequences.fasta",
            "graph with control": "../tests/lrc_kir/ctcf_with_control_sequences.fasta",
            "graph using macs reads": "../tests/lrc_kir_using_macs_reads_sequences.fasta",
            "graph using macs reads remapped": "../tests/lrc_kir/macs_reads_remapped_sequences.fasta",
        },
        "MA0139.1.meme",
        save_to_file="enrichment.png"
    )

    exit()


    plot_true_positives(
        {
            #"graph peaks2": "../tests/tmp_sequences",
            # "graph peaks on macs path": "../tests/real_data_sequences_only_macs_path.fasta",
            # "my_peaks": "../tests/ctcf_q50_with_control_sequences.fasta",
            #"my_peaks": "../tests/tmp_sequences",
            "macs reads remapped without control": "../tests/ctcf_macs_reads_remapped_without_control_sequences.fasta",
            #"macs": "CTCFpFix.fasta",
            "macs without control": "macs_no_control_peaks.fasta",
            "macs macs_using_graph_reads_peaks": "macs_using_graph_reads_peaks.fasta",
            "graph using macs reads": "../tests/linear_reads_moved_to_graph_sequences.fasta",
            "graph r0.99": "../tests/ctcf_r099_without_control_sequences.fasta",
            # "graph q30 vs q50": "../tests/ctcf_q50_vsq30_without_control_sequences.fasta",
            # "macs without control": "mac_wit.fasta",
            # "macs with control": "CTCF_filtered.fasta",
            # "graph with control": "../tests/ctcf_q50_sequences.fasta",
            "graph without control, reads outside removed": "../tests/ctcf_filtered_outside_sequences.fasta",
            "graph without control": "../tests/ctcf_q50_without_control_sequences.fasta",
            "graph without control new filtering": "../tests/ctcf_r1_sequences.fasta",
            "graph without control new filtering 2": "../tests/ctcf_filtered6_sequences.fasta",
            "graph without control subsampled": "../tests/ctcf_q60_subsampled_sequences.fasta",
            "graph with control": "../tests/ctcf_q50_with_control_2_sequences.fasta",
            #"graph using macs reads": "../tests/ctcf_macs_reads_without_control_sequences.fasta",
            #"graph using macs reads with control": "../tests/ctcf_macs_reads_with_control_sequences.fasta",
            #"graph_peaks_bugfix": "../tests/real_data_sequences_after_bugfix"
            #"macs": "CTCF_filtered.fasta"
            "macs": "CTCFpFix.fasta"
            # "new_fasta": "../tests/sequences_new_control_sample.fasta",
            # "max_instead_of_average": "../tests/real_data_sequences_max_instead_of_average.fasta",
            #"spp_from_encode": "ENCFF155DHA_filtered.fasta"

        },
        "MA0139.1.meme"
    )
