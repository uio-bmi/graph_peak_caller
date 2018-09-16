

def get_shared_ids(filename):
    return [line.split("\t")[2] for line in open(filename) if not line.startswith("#")]

if __name__ == "__main__":
    import sys
    folder = sys.argv[1]
    out_file = folder + "motif_summary.tsv"
    with open(out_file, "w") as out:
        for i in range(1, 6):
            filename = folder+ "fimo_graph_chr%s/fimo.txt" % i
            motif_ids = get_motif_ids(filename)
            out.write("%s\t%s\n" % (i, ",".join(motif_ids)))
    
    
