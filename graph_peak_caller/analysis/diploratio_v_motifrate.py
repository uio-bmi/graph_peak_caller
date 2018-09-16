import numpy as np
import matplotlib.pyplot as plt


def plot(base_name):
    def get_hist(s):
        return s["summary"][0]*s["diplo_hist"]
    motif = np.load(base_name + "/limited_summits_alignments_motif_summary.npz")
    nonmotif = np.load(base_name + "/limited_summits_alignments_nonmotif_summary.npz")
    motif_hist = get_hist(motif)
    nonmotif_hist = get_hist(nonmotif)
    cum_motif = np.cumsum(motif_hist)
    cum_nonmotif = np.cumsum(nonmotif_hist)
    cum_total = cum_motif + cum_nonmotif
    ratio = np.where(cum_total == 0, 0, cum_motif/cum_total)
    plt.plot(np.linspace(0, 1, 100), ratio, label=base_name)

if __name__ == "__main__":
    import sys
    paths = sys.argv[1].split(",")
    for path in paths:
        plot(path)
    plt.xlabel("Ratio of reads covered by diplotypes threshold")
    plt.ylabel("Motif match percentage")
    plt.legend()

    plt.show()
