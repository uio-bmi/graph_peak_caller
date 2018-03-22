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
            n_checked += 1
            true_positives.append(n_matched / n_checked)

        return true_positives


def plot_true_positives(peak_file_sets, meme_file_name, plot_title="", save_to_file=None, run_fimo=True):
    colors = ['b', 'g']
    i = 0
    font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 8}
    import matplotlib
    import numpy as np
    matplotlib.rc('font', **font)

    plt.figure(figsize=(3.0, 3.0), dpi=300)
    axis = plt.gca()
    # Major ticks every 20, minor ticks every 5
    major_ticks = np.linspace(0, 1.0, 5)
    minor_ticks = np.arange(0, 1, 20)

    possible_yticks = np.array(np.linspace(0, 1, 11))
    #ax.set_yticks(minor_ticks, minor=True)

    #ax.set_xticks(major_ticks)
    #ax.set_xticks(minor_ticks, minor=True)

    max_y = 0.0
    for name, fasta_file_name in peak_file_sets:
        matcher = MotifMatcher(fasta_file_name, meme_file_name)
        true_positives = matcher.compute_true_positives()
        n_matching = len(matcher.peaks_matching_motif)
        print("Fasta: %s" % fasta_file_name)
        n_tot = len(matcher.sorted_peaks)
        print("N tot: %d" % n_tot)
        print("True positives for %s: %d / %d = %.3f" % (name, n_matching, n_tot, n_matching/n_tot))
        plt.plot(true_positives, color=colors[i], label=name, linewidth=1.0)
        max_y = max(max_y, np.max(true_positives[2:]))
        i += 1
    #axis.set_ylim([0.75,1.0])
    axis.set_xlim(0)
    min_y = n_matching / n_tot
    idx = np.where(possible_yticks > min_y)
    yticks = possible_yticks[idx]
    yticks = yticks[np.where(yticks <= max_y)]
    axis.set_yticks(yticks)
    axis.grid(which='both')

    #plt.xlabel("Number of peaks included (of the total set of peaks sorted descending on score)")
    #plt.ylabel("Proportion of peaks enriched for motif")
    plt.title(plot_title)
    plt.legend()
    if save_to_file is not None:
        plt.savefig(save_to_file, dpi=300)
        print("Saved figure to %s" % save_to_file)
    else:
        plt.show()
