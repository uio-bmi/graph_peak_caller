class Peak(object):
    def __init__(self, chrom, start, end):
        self.chromosome = chrom
        self.start = start
        self.end = end

    @classmethod
    def from_bed_line(cls, line):
        parts = line.split()
        return cls(parts[0], int(parts[1]), int(parts[2]))


def peak_generator(bed_file_name):
    f = open(bed_file_name)
    return (Peak.from_bed_line(line) for line in f.readlines()
            if not line.startswith("#"))
