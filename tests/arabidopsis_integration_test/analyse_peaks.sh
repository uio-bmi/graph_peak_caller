

for chromosome in $(seq 1 2)
do
    graph_peak_caller get_summits -g graphs/$chromosome.nobg ${chromosome}_sequences.fasta ${chromosome}_qvalues
done

graph_peak_caller concatenate_sequence_files -s True 1,2 summits.fasta

graph_peak_caller