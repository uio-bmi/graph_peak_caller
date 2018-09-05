import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle, Polygon, FancyBboxPatch

from matplotlib.collections import PatchCollection


def square_venn(ax, gpc, macs, shared, m_gpc_u, m_gpc_s, m_macs_u, m_macs_s):
    col1 = (0.8, 0.3, 0.3)
    col2 = (0.3, 0.8, 0.3)
    total = float(gpc+macs-shared)
    r_gpc = (0, gpc, 1, col1)
    r_macs = (gpc-shared, macs, 1, col2)
    r_m_gpc_u = (0, gpc-shared, m_gpc_u/float(gpc-shared), col1)
    r_m_macs_u = (gpc, total-gpc, m_macs_u/float(macs-shared), col2)
    coords = [r_gpc, r_macs]  # , r_m_macs_u]  # , r_m_gpc_s, r_m_macs_s]
    coords = [(x/total, w/total, h, c) for x, w, h, c in coords]
    print(coords)
    pad = 0.01
    recs = [FancyBboxPatch((x+pad, 0+pad), w-2*pad, h-2*pad, alpha=0.3, facecolor=c,
                           linestyle="dashed", boxstyle="round,pad=0.01")
            for x, w, h, c in coords]
    # recs = [FancyBboxPatch((x, 0), w, h, alpha=0.3, facecolor=c,
    #                        linestyle="dashed", boxstyle="round, pad=0.0")
    #         for x, w, h, c in coords]

    p_gpc_u = m_gpc_u/float(gpc-shared)
    p_gpc_s = m_gpc_s/float(shared)
    p_macs_u = m_macs_u/float(macs-shared)
    p_macs_s = m_macs_s/float(shared)
    xy_gpc_s = [(0, 0),
                (gpc/total, 0),
                (gpc/total, p_gpc_s),
                ((gpc-shared)/total, p_gpc_s),
                ((gpc-shared)/total, p_gpc_u),
                (0, p_gpc_u)]
    xy_macs_s = [(1, 0),
                 (1-macs/total, 0),
                 (1-macs/total, p_macs_s),
                 (gpc/total, p_macs_s),
                 (gpc/total, p_macs_u),
                 (1, p_macs_u)]
    polygs = []
    polygs.append(
        Polygon(xy=xy_gpc_s, facecolor=col1, alpha=0.7))
    polygs.append(
        Polygon(xy=xy_macs_s, facecolor=col2, alpha=0.7))

    # r_m_gpc_s = ((gpc-shared, 0), (gpc, u_gpc_s/float(shared)), 
    # r_m_gpc_s = (gpc-shared, shared, m_gpc_s/float(shared), "g")
    # r_m_macs_s = (gpc-shared, shared, m_macs_s/float(shared), "r")
    # collection = PatchCollection(recs, alpha=0.3)
    # collection.set_array(np.array([0, 0.25, 0.5, 0.75, 0.5, 0.75]))
    # collection.set_array(np.linspace(0, 1, len(recs)))
    for rec in recs + polygs:
        ax.add_patch(rec)

    return (gpc-shared)/(2*total), (gpc-0.5*shared)/total, (gpc+(macs-shared)/2.)/total
    # ax.add_collection(collection)

def column_venn(ax, gpc, macs, shared, m_gpc_u, m_gpc_s, m_macs_u, m_macs_s, both_motif):
    total = gpc+macs-shared
    xs = [0, (gpc-shared)/total, (gpc+macs-2*shared)/total]
    widths = [(gpc-shared)/total, (macs-shared)/total, shared/total]
    print(xs, widths)
    heights = [m_gpc_u/(gpc-shared), m_macs_u/(macs-shared), both_motif/shared]
    colors = ["b", "#db9200", "g"]
    big_rectangles = [Rectangle((x, 0), w, 1, facecolor=c, alpha=0.3)
                      for x, w, c in zip(xs, widths, colors)]
    small_rectangles = [Rectangle((x,0), w, h, facecolor=c, alpha=0.3)
                        for x, w, h, c in zip(xs, widths, heights, colors)]
    last_h = (m_gpc_s+m_macs_s-2*both_motif)/shared
    ws = [widths[2]*m_gpc_s/(m_gpc_s+m_macs_s), widths[2]*m_macs_s/(m_gpc_s+m_macs_s)]
    last_recs = [Rectangle((xs[2], heights[2]), ws[0], last_h,
                           facecolor="b", alpha=0.6),
                 Rectangle((xs[2]+ws[0], heights[2]), ws[1], last_h,
                           facecolor="#e0b721", alpha=0.6)]

    for rec in big_rectangles+small_rectangles+last_recs:
        ax.add_patch(rec)
    return [x+w*0.5 for x, w in zip(xs, widths)]


def from_csv(filename):
    lines = [line.split() for line in open(filename)]
    counts = [float(n) for n in lines[1][1:]]
    print(counts)
    names = lines[0][1:]
    d = dict(zip(names, counts))
    numbers = counts[:3]+[
        d["MOTIF_UNIQUE_GPC"],
        d["MOTIF_SHARED_GPC"],
        d["MOTIF_UNIQUE_MACS"],
        d["MOTIF_SHARED_MACS"],
        d["MOTIF_BOTH"]]
    return numbers


def save_venn(numbers, out_name):
    fig, ax = plt.subplots(1)
    xs = column_venn(ax, *numbers)
    #plt.text(xs[0], 0.8, "GPC", ha="center")
    #plt.text(xs[1], 0.8, "MACS", ha="center")
    #plt.text(xs[2], 0.8, "SHARED", ha="center")
    try:
        # Try hacky way to get title
        title = out_name.split("/")[-1].split(".")[0].split("_")[1]
    except IndexError:
        title = "Untitled %s" % out_name

    plt.title(title)
    plt.savefig(out_name+"_venn.pdf", bbox_inches='tight')


def save_venn_from_csv(filename, out_name):
    numbers = from_csv(filename)
    save_venn(numbers, out_name)

def save_venn_from_csvs(filenames, out_name):
    #fig, axes = plt.subplots(2,3)
    #axes = [ax for row in axes for ax in row]
    numbers_list = [from_csv(filename) for filename in filenames]
    #for numbers, ax in zip(numbers_list, axes):
    #    column_venn(ax, *numbers)
    total = [sum(elems) for elems in zip(*numbers_list)]
    #print(total)
    #column_venn(axes[-1], *total)
    font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 30}
    import matplotlib
    matplotlib.rc('font', **font)
    plt.figure(figsize=(15, 15.0), dpi=300)
    ax = plt.gca()
    column_venn(ax, *total)
    plt.show()
    plt.savefig(out_name)
    # plt.text(xs[0], 0.8, "GPC", ha="center")
    # plt.text(xs[1], 0.8, "MACS", ha="center")
    # plt.text(xs[2], 0.8, "SHARED", ha="center")
    # plt.savefig(out_name+"_venn.pdf", bbox_inches='tight')
    # numbers = from_csv(filename)
    # save_venn(numbers, out_name)

if __name__ == "__main__":
    import sys
    filenames = sys.argv[1]
    save_venn_from_csvs(filenames.split(","), sys.argv[2])
