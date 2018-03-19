graph_dir="graphs/"
motif="motif.meme"

# Run fimo for each chromosome
for chromosome in $(seq 1 2)
do
    fimo -oc fimo_$chromosome $motif ${chromosome}_sequences.fasta
done

for chromosome in $(seq 1 2)
do
    graph_peak_caller diffexpr -g $graph_dir/$chromosome.nobg $chromosome fimo_$chromosome/fimo.txt
    cat ${chromosome}_diffexpr.fasta >> diff_expr.fasta
done

for chromosome in $(seq 1 2)
do
    graph_peak_caller get_summits -g graphs/$chromosome.nobg ${chromosome}_sequences.fasta ${chromosome}_qvalues
done

#graph_peak_caller analyse_peaks_whole_genome 1,2 ./ graphs/ results
