from offsetbasedgraph.tracevariants import *
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.stats import binom_test
import motif_models

count_mismatch = lambda result: result.total - (result.A_count + result.B_count)

def check_is_diplo(result):
    return result.total == result.A_count+result.B_count


def _get_binomial_ps(mismatches, total, rate):
    return np.array([binom_test(missees, n, p=rate, alternative="greater")
                     for missees, n in zip(mismatches, total)])


def get_motif_plot(motif_results, all_results):
    motif_mismatches = np.array([count_mismatch(result) for result in motif_results])
    all_mismatches = np.array([count_mismatch(result) for result in all_results])
    motif_total = np.array([result.total for result in motif_results])
    all_total = np.array([result.total for result in all_results])
    mismatch_rate = np.median(all_mismatches/all_total)
    threshold = 1
    motif_ps = np.where(motif_total>0, 1-motif_mismatches/motif_total, 0)
    all_ps = np.where(all_total>0, 1-all_mismatches/all_total, 0)
    motif_counts = np.cumsum(np.histogram(motif_ps[::-1], 100, (0, 1))[0])
    all_counts = np.cumsum(np.histogram(all_ps[::-1], 100, (0, 1))[0])
    sorted_all = np.sort(1-all_ps)
    sorted_motif = np.sort(1-motif_ps)
    motif_plot = np.digitize(sorted_all, sorted_motif)/(1+np.arange(sorted_all.size))
    idxs = np.r_[np.flatnonzero(np.diff(sorted_all)), motif_plot.size-1]
    motif_plot = motif_plot[idxs]
    return idxs, motif_plot
    # return sorted_all[idxs], motif_plot


def compare_two(motif_results_A, all_results_A, motif_results_B, all_results_B):
    motif_plot_A = get_motif_plot(motif_results_A, all_results_A)
    motif_plot_B = get_motif_plot(motif_results_B, all_results_B)
    # args = np.argsort(np.r_[motif_plot_A[0], motif_plot_B[0]])
    # B_idxs = args >= motif_plot_A[0].size
    # A_idxs = ~B_idxs
    # B_idxs = np.flatnonzero(B_idxs)
    # A_idxs = np.flatnonzero(A_idxs)
    A_idxs = motif_plot_A[0]
    B_idxs = motif_plot_B[0]
    print(A_idxs.shape, motif_plot_A[1].shape)
    print(B_idxs.shape, motif_plot_B[1].shape)
    plt.plot(A_idxs, motif_plot_A[1], label="graph")
    plt.plot(B_idxs, motif_plot_B[1], label="macs")
    plt.legend()
    # plt.show()
    # plt.plot(motif_plot_A[1], label="graph")
    # plt.plot(motif_plot_B[1], label="macs")

def analyze_result(motif_results, all_results, col="b"):
    motif_mismatches = np.array([count_mismatch(result) for result in motif_results])
    all_mismatches = np.array([count_mismatch(result) for result in all_results])
    motif_total = np.array([result.total for result in motif_results])
    all_total = np.array([result.total for result in all_results])
    mismatch_rate = np.median(all_mismatches/all_total)
    threshold = 1
    print(mismatch_rate)
    # mismatch_rate = np.sum(all_mismatches)/np.sum(all_total)
    # print(mismatch_rate)
    motif_ps = np.where(motif_total>0, 1-motif_mismatches/motif_total, 0)
    all_ps = np.where(all_total>0, 1-all_mismatches/all_total, 0)
    # motif_ps = _get_binomial_ps(motif_mismatches, motif_total, mismatch_rate)
    # all_ps = _get_binomial_ps(all_mismatches, all_total, mismatch_rate)
    # plt.hist(_motif_ps, 20, (0, 1))
    # plt.show()
    # plt.hist(_all_ps, 20, (0, 1))
    # plt.show()
    motif_counts = np.cumsum(np.histogram(motif_ps[::-1], 100, (0, 1))[0])
    # plt.hist(motif_ps, 200)
    # plt.show()
    all_counts = np.cumsum(np.histogram(all_ps[::-1], 100, (0, 1))[0])
    # plt.hist(all_ps, 200)
    # plt.show()
    sorted_all = np.sort(1-all_ps)
    sorted_motif = np.sort(1-motif_ps)
    motif_plot = np.digitize(sorted_all, sorted_motif)/(1+np.arange(sorted_all.size))
    idxs = np.r_[np.flatnonzero(np.diff(sorted_all)), motif_plot.size-1]
    # split = np.flatnonzero(sorted_all[idxs] > (1-threshold))[0]
    motif_plot = motif_plot[idxs]
    # plt.plot(idxs[:split], motif_plot[:split], col)
    # plt.plot(idxs[split:], motif_plot[split:], col)
    plt.plot(idxs, motif_plot, col)
    # plt.plot(idxs, sorted_all[idxs]*np.max(motif_plot), "g")
    # plt.plot(sorted_all*np.max(motif_plot)/np.max(sorted_all))
    a=np.count_nonzero(motif_mismatches==0)
    b=np.count_nonzero(all_mismatches==0)
    print("#PRE:",motif_ps.size, all_ps.size, motif_ps.size/all_ps.size)
    print("#POST:", a, b, a/float(b))
    # print("#POST:", (motif_ps>=threshold).sum(), (all_ps>=threshold).sum(), (motif_ps>=threshold).sum()/float((all_ps>=threshold).sum()))
    return (a, b)

def count_diplo(results):
    return sum(check_is_diplo(result) for result in results)


def parse_results(filename):
    return [(line.split("\t")[0], eval(line.split("\t")[1]))
            for line in open(filename) if not line.startswith("#")]

def run_comparison(folder, unique=True):
    print(folder)
    if unique:
        macs_results = {i: parse_results("%s%s_macs_unique_summits_alignments_diplotypes.tsv" % (folder, i)) for i in range(1, 6)}
        is_unique = lambda peak, i: peak in graph_unique[str(i)]
    else:
        macs_results = {i: parse_results("%s%s_macs_all_summits_alignments_diplotypes.tsv" % (folder, i)) for i in range(1, 6)}        
        is_unique= lambda peak, i: True
    graph_results = {i: parse_results("%s%s_limited_summits_alignments_diplotypes.tsv" % (folder, i)) for i in range(1, 6)}
    macs_motifs = {line.split("\t")[0]: {peak.strip() for peak in line.split("\t")[1].split(",")} for line in open(folder+"/fimo_macsmotif_summary.tsv")}
    graph_motifs = {line.split("\t")[0]: {peak.strip() for peak in line.split("\t")[1].split(",")} for line in open(folder+"/fimo_graphmotif_summary.tsv")}
    graph_unique = {line.split("\t")[0]: {peak.strip() for peak in line.split("\t")[1].split(",")} for line in open(folder+"/unqiue_summary.tsv")}

    macs_non_motif_results = [result for i, results in macs_results.items()
                              for peak, result in results if  peak not in macs_motifs[str(i)]]
    macs_motif_results = [result for i, results in macs_results.items()
                     for peak, result in results if peak in macs_motifs[str(i)]]
    graph_non_motif_results = [result for i, results in graph_results.items()
                         for peak, result in results if is_unique(peak, i) and peak not in graph_motifs[str(i)]]
    graph_motif_results = [result for i, results in graph_results.items()
                           for peak, result in results if is_unique(peak, i) and peak in graph_motifs[str(i)]]
    # compare_two(graph_motif_results, graph_motif_results+graph_non_motif_results,
    # macs_motif_results, macs_motif_results+macs_non_motif_results)
    # plt.show()
    return graph_motif_results, graph_non_motif_results, macs_motif_results, macs_non_motif_results
    motif_models.simple_model(graph_motif_results, graph_non_motif_results,
                              macs_motif_results, macs_non_motif_results)


def compare_folders(folders):
    results = [[], [], [], []]
    counts = []
    tots = []
    for folder in folders:
        res = run_comparison(folder)

        for i in range(4):
            results[i].extend(res[i])
        print([len(r) for r in res])
        print(folder)
        counts.append([count_diplo(r) for r in res])
        tots.append([len(r) for r in res])
        # compare_two(res[0], res[0]+res[1], res[2], res[2]+res[3])
        # plt.show()
    print(counts)
    print(tots)
    tots = np.array(tots)
    counts = np.array(counts)
    inv = tots-counts
    props1 = counts[:, 0]/(counts[:, 0]+counts[:, 1])
    props2 = counts[:, 2]/(counts[:, 2]+counts[:, 3])
    props3 = inv[:, 0]/(inv[:, 0]+inv[:, 1])
    props4 = inv[:, 2]/(inv[:, 2]+inv[:, 3])
    print("graph1", props1)
    print("graph2", props3)
    print("graphdiff", props1-props3)
    print("macs1", props2)
    print("macs2", props4)
    print("macsdiff", props2-props4)

    # props = counts/tots
    # print(tots.T)
    # print(counts.T)
    # print(props.T)
    res = results
    compare_two(res[0], res[0]+res[1], res[2], res[2]+res[3])
    plt.show()
    motif_models.simple_model(*results)
        

def run(folder, includes="all"):
    #interval_name = "macs_unique_graph_"
    interval_name= ""
    rs = {i: parse_results("%s%s_%slimited_summits_alignments_diplotypes.tsv" % (folder, i, interval_name)) for i in range(1, 6)}
    with_macs = includes == "macs"
    if includes == "macs":
        summary_name = "motif_summary_graph_matching_macs.tsv"
        includes = "shared"
        macs_motifs = {line.split("\t")[0]: {peak.strip() for peak in line.split("\t")[1].split(",")} for line in open(folder+"/" + summary_name)}        
    summary_name = "fimo_graphmotif_summary.tsv"
    motifs = {line.split("\t")[0]: {peak.strip() for peak in line.split("\t")[1].split(",")} for line in open(folder+"/" + summary_name)}        
    print(sum(len(m) for m in motifs.values()))
    # summary_name = "fimo_macs_uniquemotif_summary.tsv"

    unique = {line.split("\t")[0]: {peak.strip() for peak in line.split("\t")[1].split(",")} for line in open(folder+"/unqiue_summary.tsv")}

    # for chrom, motif in motifs.items():
    #     in_chrom = {result[0] for result in rs[int(chrom)]}
    #     missing = [m for m in motif if m not in in_chrom]
    #     assert not missing, missing
    #     # assert all(m in in_chrom for m in motif)
    # 
    # for chrom, _unique in unique.items():
    #     in_chrom = {result[0] for result in rs[int(chrom)]}
    #     assert all(u in in_chrom for u in _unique)

    filter_funcs = {"all": lambda peak, i: True,
                    "unique": lambda peak, i: peak in unique[str(i)],
                    "shared": lambda peak, i: peak not in unique[str(i)]}
    
    non_motif_results = [result for i, results in rs.items()
                         for peak, result in results if filter_funcs[includes](peak, i) and peak not in motifs[str(i)]]
    motif_results = [result for i, results in rs.items()
                     for peak, result in results if filter_funcs[includes](peak, i) and peak in motifs[str(i)]]
    all_results = motif_results + non_motif_results
    plt.title(folder)
    print(folder)
    res = analyze_result(motif_results, all_results)
    if with_macs:
        macs_motif_results = [result for i, results in rs.items()
                         for peak, result in results if filter_funcs[includes](peak, i) and peak in macs_motifs[str(i)]]

        analyze_result(macs_motif_results, all_results, "r")
    # plt.show()
    return res


def run_func(folder, interval_name, func):
    includes = "all"
    rs = {i: parse_results("%s%s_%s_diplotypes.tsv" % (folder, i, interval_name)) for i in range(1, 6)}
    summary_name = "fimo_graphmotif_summary.tsv"
    motifs = {line.split("\t")[0]: {peak.strip() for peak in line.split("\t")[1].split(",")} for line in open(folder+"/" + summary_name)}        
    unique = {line.split("\t")[0]: {peak.strip() for peak in line.split("\t")[1].split(",")} for line in open(folder+"/unqiue_summary.tsv")}
    shared_results = [result for i, results in rs.items()
                      for peak, result in results if peak not in unique[str(i)]]

    unique_results = [result for i, results in rs.items()
                     for peak, result in results if peak in unique[str(i)]]
    print(folder)
    print("Unique")
    func(unique_results)
    print("Shared")
    func(shared_results)
    print("All")
    func(unique_results+shared_results)


def count_haplo(results):
    haplo = sum(result.A_count == result.total for result in results)
    t = len(results)
    if t == 0:
        print(haplo, t, "NA")
    else:
        print(haplo, t, haplo/t)

if __name__ == "__main__":
    if sys.argv[1] == "haplo":
        [run_func(folder+"/", sys.argv[2], count_haplo) for folder in sys.argv[3:]]
        exit()
    if sys.argv[1] == "compare":
        compare_folders(sys.argv[2:])
        # [run_comparisons(folder+"/") for folder in sys.argv[2:]]
        exit()



    folders = sys.argv[1:]
    output = [run(folder+"/", includes="unique") for folder in folders]
    s = np.sum(output, axis=0)
    results = np.array(output).T
    print(" ".join(folder.split("/")[1].split("_")[1] for folder in folders))
    for row in results:
        print(" ".join(str(n) for n in row))
    print(" ".join(str(n) for n in results[0]/results[1]))
    print(s[0], s[1], s[0]/s[1])
    plt.show()
