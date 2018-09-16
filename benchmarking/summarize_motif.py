

def get_motif_ids(filename):
    return {line.split("\t")[2] for line in open(filename) if not line.startswith("#")}

def run(folder, fimo_folder):
    out_file = folder + fimo_folder + "motif_summary.tsv"
    counter = 0
    with open(out_file, "w") as out:
        for i in range(1, 6):
            filename = folder+ fimo_folder + "_chr%s/fimo.txt" % i
            motif_ids = get_motif_ids(filename)
            counter += len(motif_ids)
            print(i, len(motif_ids))
            out.write("%s\t%s\n" % (i, ",".join(motif_ids)))
    print("Number fo motifs %s " % counter)
    

if __name__ == "__main__":
    import sys
    fimo_folder = sys.argv[1]
    folders = sys.argv[2:]
    [run(folder, fimo_folder) for folder in folders]
    
