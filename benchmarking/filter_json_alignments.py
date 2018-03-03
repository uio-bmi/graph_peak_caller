import re
import logging
import sys

ids_to_remove = set()
f = open(sys.argv[1])
i = 0
for line in f:
    if i % 500000 == 0:
        print("%d ids read" % i)
    i += 1
    ids_to_remove.add(line.rstrip())

f.close()
print("Done  reading ids. Now filtering")

outfile = open(sys.argv[3], "w")
f = open(sys.argv[2])
regex = re.compile(r"name\": \"(.*?)\"")
i = 0
n_removed = 0
for line in f:
    if i % 500000 == 0:
        print("%d reads processed. %d removed so far." % (i, n_removed))
    i += 1
    groups = regex.findall(line)
    name = None
    for n in groups:
        if len(n) > 3:
            name = n.split()[0]

    if name is None:
        print("No name")
        continue

    #name = groups[0].split()[0]
    name = name.replace("/1", "")
    if name not in ids_to_remove:
        outfile.writelines([line])
    else:
        n_removed += 1
outfile.close()
f.close()
print("In total %d reads removed." % n_removed)
print("Resulting reads written to %s" % sys.argv[3])

