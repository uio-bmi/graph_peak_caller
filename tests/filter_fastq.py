def remove_reads_in_other_fastq(input_fastq, remove_fastq, out_file_name):

    out = open(out_file_name, "w")
    # Get seq ids to remove
    remove_ids = set()
    i = 0
    for line in open(remove_fastq):
        if i % 500000 == 0:
            print(i)
        i += 1
        if line.startswith("@"):
            remove_ids.add(line)

    i = 0
    write_line = True
    for line in open(input_fastq):
        if i % 100000 == 0:
            print(i)
        i += 1
        line = line.replace("/1", "")
        if line.startswith("@"):
            if line in remove_ids:
                write_line = False
            else:
                write_line = True

        if write_line:
            out.writelines(["%s" % line])


if __name__ == "__main__":
    remove_reads_in_other_fastq("orig.fastq", "filtered_outside_mhc.fastq", "final.fastq")




