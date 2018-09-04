import json
import pandas
import numpy as np


def peaks_from_fasta(filename):
    with open(filename) as f:
        parts = (line.split(maxsplit=1) for line in f if line.startswith(">"))
        entries = [json.loads(part[1]) for part in parts]
        return pandas.DataFrame(entries)


def df_from_fimo(filename):
    return pandas.read_csv(filename, sep="\t").set_index("sequence_name")

def df_from_subgraphs(filename):
    count_edges = lambda x: np.count_nonzero(x.data)
    count_nodes = lambda x: x.shape[0]-1
    subgraphs = {name: graph[()] for name, graph in np.load(filename).items()}
    return pandas.DataFrame([{"name": name, "edges": count_edges(graph), "nodes": count_nodes(graph)}
                             for name, graph in subgraphs.items()]).set_index("name")
    
def df_from_setsummary(filename):
    names, parts = zip(*(line.split(", (", 1) for line in open(filename) if not line.startswith("#")))
    # names = ["peak"+name for name in names]
    parts = [part.split(", ", 2) for part in parts]
    ratios = [(float(a), float(b)) for a, b, _ in parts]
    df = pandas.DataFrame([{"coverage": r[0], "reads": r[1]}
                           for r in ratios])
    return df

def combine_dfs(folder, name):
    peaks = peaks_from_fasta(folder + name + "_sequences_summits_unique.fasta")
    fimo = df_from_fimo(folder + "fimo_graph_chr"+ name + "/fimo.txt")
    tmp = df_from_subgraphs(folder+name+"_sub_graphs.graphs.npz")
    coverage = df_from_setsummary(folder + name + "_sequences_summits_unique.setsummary")
    combined = pandas.concat([peaks, coverage], axis=1, join="inner")
    combined = combined.set_index("unique_id")
    subgraphs = pandas.concat([combined, tmp], axis=1, join="inner")
    match_ids = set(fimo.index.unique())
    subgraphs["MATCH"] = [name in match_ids for name in subgraphs.index.values]
    print(subgraphs.shape)
    return subgraphs
    matches = subgraphs.loc[subgraphs["MATCH"]]
    nonmatches = subgraphs.loc[~subgraphs["MATCH"]]
    print("MATCHES", (matches["edges"]-matches["nodes"]).mean(), matches.shape[0])
    print("NONMAtCHES", (nonmatches["edges"]-nonmatches["nodes"]).mean(), nonmatches.shape[0])


def motif_haplotype_analysis(folder, name):
    coverage = df_from_setsummary(folder + name + "_sequences_summits_unique.setsummary")


if __name__ == "__main__":
    dfs = [combine_dfs("./", str(c)) for c in range(1, 6)]
    combined = pandas.concat(dfs, axis=0, keys=list(range(1, 6)))
    print(combined["MATCH"].mean())
    print(combined["MATCH"].sum())
    tmp = combined.loc[(combined["coverage"]/combined["reads"])>=0.8]
    print(tmp)
    print(tmp["MATCH"].mean())
    print(combined.loc[(combined["coverage"]/combined["reads"])<0.8]["MATCH"].mean())
    # for c in (str(i) for i in range(1,6)):
    #     combine_dfs("./", # "/home/ivar/dev/graph_peak_caller/benchmarking/data/ARABIDOPSIS_ERF115_SRR931836/1/",
    #                 c)

    # peaks = peaks_from_fasta("/home/ivar/dev/graph_peak_caller/benchmarking/data/ARABIDOPSIS_ERF115_SRR931836/1/unique_graph.fasta")
    # fimo = df_from_fimo("/home/ivar/dev/graph_peak_caller/benchmarking/data/ARABIDOPSIS_ERF115_SRR931836/1/fimo_unique_graph/fimo.txt")
    # subgraphs = df_from_subgraphs("/home/ivar/dev/graph_peak_caller/benchmarking/data/ARABIDOPSIS_ERF115_SRR931836/1/1_sub_graphs.graphs.npz")
    # print(subgraphs)

