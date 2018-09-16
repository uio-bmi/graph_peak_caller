null = None

def get_unique_ids(filename):
    return [eval(line.strip())["unique_id"] for line in open(filename)]

def run(folder):
    out_file = folder + "unqiue_summary.tsv"
    unique_count = 0
    with open(out_file, "w") as out:
        for i in range(1, 6):
            filename = folder+ "/%s_sequences_limited_summits_unique.intervalcollection" % i
            unique_ids = get_unique_ids(filename)
            unique_count+=len(unique_ids)
            out.write("%s\t%s\n" % (i, ",".join(unique_ids)))

    print("Unique Ids", unique_count)

if __name__ == "__main__":
    import sys
    folders = sys.argv[1:]
    [run(folder) for folder in folders]
