import sys
from graph_peak_caller.peakcollection import PeakCollection


def parse_set_file(file_name):
    lines = (line for line in open(file_name) if not line.startswith("#"))
    parts = (line.split(", (", 1) for line in lines)
    parts = (part[1].split(", ", 2) for part in parts)
    return  ((int(part[0]), int(part[1]))  for part in parts)

def parse_dip_set_file(file_name):
    lines = (line for line in open(file_name) if not line.startswith("#"))
    parts = (line.split(", (", 1) for line in lines)
    parts = (part[1].split(", ", 3) for part in parts)
    return  ((int(part[0]), int(part[1]), int(part[2]))  for part in parts)

def check_unique(peak_file_name, motif_file_name, set_file_name, filter_unique):
    is_dip = "_dip." in set_file_name
    peaks = PeakCollection.from_file(peak_file_name, True)
    ids = {peak.unique_id for peak in peaks}
    motifs = PeakCollection.from_file(motif_file_name, True)
    unique_motifs = {i for i, motif in enumerate(motifs) if motif.unique_id in ids}
    parser = parse_dip_set_file if is_dip else parse_set_file
    counts = [count for i, count in enumerate(parser(set_file_name)) if (not filter_unique) or (i in unique_motifs)]
    N = len(counts)
    print(N, motif_file_name)
    if is_dip:
        successes = len([count for count in counts if (count[0]+count[1]) >= count[2]])
    else:
        successes = [count for count in counts if count[0] >= count[1]]
        print(successes)
        successes = len(successes)
    return successes, N


if __name__ == "__main__":
    folder, name, chromosomes, filter_unique, strictness = sys.argv[1:6]
    set_filenames = (folder+chromosome+"_" + name + strictness + ".setsummary" for chromosome in chromosomes.strip().split(","))
    motif_filenames = (folder+chromosome+"_" + name + ".intervalcollection" for chromosome in chromosomes.strip().split(","))
    unique_filenames = (folder+chromosome+"_" + "sequences_summits_unique.intervalcollection" for chromosome in chromosomes.strip().split(","))
    total, s = (0, 0)
    for peak_filename, motif_filename, set_filename in zip(unique_filenames, motif_filenames, set_filenames):
        success, N = check_unique(peak_filename, motif_filename, set_filename, filter_unique=int(filter_unique))
        total += N
        s += success
    print("#", total, s, s/total)
    # 
    # 
    # firstlines = (open(filename).readline() for filename in filenames)
    # successes, totals = zip(*(line.split()[1].strip().split("/") for line in firstlines))
    # print(sum(int(s) for s in successes)/sum(int(t) for t in totals))


