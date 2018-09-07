import sys
import logging

bed = open(sys.argv[1])
found = set()

for line in bed:
    l = line.split() 
    chrom = l[0]
    start = l[1]
    strand = l[5]
    if int(l[2]) <= int(l[1]):
        logging.warning("End %d <= start %s. Skipping" % (int(l[2]), start))
        continue
    pos = chrom + start + strand
    if pos in found:
        continue

    found.add(pos)
    print(line.strip())
      
