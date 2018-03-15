graph_dir="graphs/"

for chromosome in $(seq 1 2)
do

    graph_peak_caller callpeaks_whole_genome $chromosome \
        -d $graph_dir -s arabidopsis_sample/sample_ -n "" -f 230 -r 50 \
        -p True \
        -u 250000 \
        -g 135000000
done
wait

for chromosome in $(seq 1 2)
do
    graph_peak_caller callpeaks_whole_genome_from_p_values $chromosome \
            -d $graph_dir -f 230 -r 50

done

wait

# Run fimo for each chromosome
for chromosome in $(seq 1 2)
do
    fimo -oc fimo_$chromosome $motif ${chromosome}_sequences.fasta > fimo_log.txt 2>&1
done

for chromosome in $(seq 1 2)
do
    graph_peak_caller diffexpr -g $graph_dir/$chromosome.nobg $chromosome fimo_$chromosome/fimo.txt $graph_dir/$chromosome.json
    cat ${chromosome}_diffexpr.fasta >> diff_expr.fasta
done
