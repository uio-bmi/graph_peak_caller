import sys

bed = open(sys.argv[1])
found = set()

for line in bed:
    l = line.split() 
    chrom = l[0]
    start = l[1]
    strand = l[5]
    pos = chrom + start + strand
    if pos in found:
        continue

    found.add(pos)
    print(line.strip())
      
