import sys

narrowpeak_file = open(sys.argv[1])
n_bp_each_side = int(sys.argv[2])
out_file = open(sys.argv[3], "w")

for line in narrowpeak_file:
    if line.startswith("track"):
        out_file.writelines([line])
        continue
    l = line.split()
    chrom = l[0]
    start = int(l[1])
    end = int(l[2])
    peak = int(l[9])

    new_start = max(start, start + peak - n_bp_each_side)
    new_end = min(end, start + peak + n_bp_each_side)

    assert new_end - new_start <= 2*n_bp_each_side + 1

    out_file.writelines(["chr%s\t%d\t%d\t%s\t%s\t.\t%s\t%s\t%s\t%s\n" % (chrom, new_start, new_end, l[3], l[4], l[6], l[7], l[8], (start + peak) - new_start)])


out_file.close()
narrowpeak_file.close()

