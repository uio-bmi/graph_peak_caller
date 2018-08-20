import re
import sys
import numpy as np

def get_alignment_ids(data_dir, chromosomes):
    alignments = set()
    for chromosome in chromosomes:
        file_name = data_dir + "/" + "peaks1_alignments_chr" + chromosome + ".txt"
        with open(file_name) as f:
            for line in f:
                alignments.add(line.strip())

    print("Found %d alignments" % len(alignments))
    return alignments


def get_linear_alignments(data_dir, alignments, linear_alignments_out_file_name="out.tmp"):
    scores = []
    optimal_hits = []
    linear_alignments = {}
    
    out_file = open(linear_alignments_out_file_name, "w")

    with open(data_dir + "/alignments.sam") as f:
        i = 0
        for line in f:
            if line.startswith("@"):
                print("Skipping header")
                continue
            l = line.split()
            if l[0] in alignments:
                scores.append(int(l[4]))
                if len(l) >= 14:
                    optimal_hits.append(int(l[13].replace("X0:i:", "")))
                    linear_alignments[l[0]] = (int(l[2]), int(l[3]))
                    
                    out_file.writelines(line)

            if i % 1000000 == 0:
                print("Line %d" % i) 
            i+= 1

    print("Number of linear alignments: %d" % len(scores))
    print("Number of linear alignments that aligned: %d" % len(linear_alignments))
    print(np.mean(scores))
    print(np.mean(optimal_hits))
    out_file.close()
    return linear_alignments


def get_graph_alignments(data_dir, alignment_ids):
   
    graph_alignments = {}
    regex = re.compile(r"name\":\"(.*?)\"")
    regex_refpos = re.compile(r"\"refpos\":\[\{\"offset\":\"([0-9]*)\",\"name\":\"([0-9])\"\}\]")
    regex_reverse_refpos = re.compile(r"\"refpos\":\[\{\"offset\":\"([0-9]*)\",\"is_reverse\":true,\"name\":\"([0-9])\"\}\]")
    with open(data_dir + "/filtered_low_qual_reads_removed.json") as f:
        i = 0
        for line in f:
            groups = regex.findall(line) 
            name = None
            for n in groups:
                if len(n) > 3:
                    name = n.split()[0].strip()
            
            if i % 1000000 == 0:
                print("Read %d graph alignments" % i)
            i += 1

            if name not in alignment_ids:
                continue            
 
            refpos_groups = regex_refpos.findall(line)
            reverse_refpos_groups = regex_refpos.findall(line)
            
            if len(refpos_groups) == 0:
                assert len(reverse_refpos_groups) > 0
                refpos_groups = reverse_refpos_groups 
             
            if name in alignment_ids:
                graph_alignments[name] = (int(refpos_groups[0][1]), int(refpos_groups[0][0]))
            

    print("Found %d graph alignments" % len(graph_alignments))
    return graph_alignments

def compare_graph_and_linear(linear_alignments, graph_alignments):
    
    not_in_linear = 0
    different_chrom = 0
    different_pos = 0
    same_pos = 0

    for name in graph_alignments.keys():

        if name not in linear_alignments:
            not_in_linear += 1 
            continue
        
        linear_chrom = linear_alignments[name][0]
        graph_chrom = graph_alignments[name][0]
  
        if linear_chrom != graph_chrom:
            print("Linear chrom %s != graph chrom %s" % (linear_chrom, graph_chrom))
            different_chrom += 1
            continue

        linear_pos = linear_alignments[name][1]
        graph_pos = graph_alignments[name][1]
        if abs(linear_pos - graph_pos) > 100:
            different_pos += 1
        else:
            same_pos += 1
            #print("Same: %s --- %s" % (linear_alignments[name], graph_alignments[name]))
    
    print("Not in linear: %d" % not_in_linear)
    print("Different chrom: %d" % different_chrom)
    print("Different pos: %d" % different_pos)
    print("Same position: %d" % same_pos)



if __name__ == "__main__":
    data_dir = sys.argv[1]
    chromosomes = sys.argv[2].split(",") 
    alignments = get_alignment_ids(data_dir, chromosomes)
    linear_alignments = get_linear_alignments(data_dir, alignments)
    graph_alignments = get_graph_alignments(data_dir, alignments)
    compare_graph_and_linear(linear_alignments, graph_alignments)


