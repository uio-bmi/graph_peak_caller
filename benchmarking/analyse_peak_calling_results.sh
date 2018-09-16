#!/usr/bin/env bash

encode_id=$1
replicate=$2
tf=$3
chromosomes=$4
motif_url=$5
data_dir=$6
linear_genome_fasta_file=$7
genome_size=$8
peaks_file_name=$9
summit_window_size=${10}

base_dir=$(pwd)
n_threads=$(grep -c ^processor /proc/cpuinfo)

echo "Running analyse_peak_calling_results.sh"
work_dir=data/${tf}_${encode_id}/${replicate}
cd $work_dir

echo "Changed dir to $work_dir"

echo "Will use window size $summit_window_size"

background_model_file="$data_dir/background.model"

echo "Will use fimo background model file from $background_model_file"

# fimo 


# Merge all graph peak caller peaks into one file in order to 
# sort on score and choose the same amount of peaks as max
sort -n -r -k 9 $peaks_file_name > macs_peaks_sorted.narrowPeak
graph_peak_caller concatenate_sequence_files -f [chrom]_sequences.fasta $chromosomes all_graph_peaks.fasta
n_macs=$(wc -l $peaks_file_name | cut -d " " -f 1)
n_graph=$(wc -l all_graph_peaks.fasta | cut -d " " -f 1)
n_graph=$(($n_graph/2))

if (( n_macs > n_graph ));
 then
    limit=$n_graph
else
    limit=$n_macs
fi

limit="100000000"  # No limiting
fasta_limit=$(($limit*2))
echo "Limiting all peaks to $limit (fasta limit $fasta_limit)"
head -n $limit all_graph_peaks.intervalcollection > all_graph_peaks_limited.intervalcollection
head -n $limit macs_peaks_sorted.narrowPeak > all_macs_peaks_limited.narrowPeak

# Split graph peaks by chromosome again
graph_peak_caller split_peaks_by_chromosome all_graph_peaks_limited.intervalcollection $chromosomes limited.intervalcollection

# Get sequences for each graph chromosome
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    graph_peak_caller peaks_to_fasta $data_dir/$chromosome.nobg.sequences ${chromosome}_limited.intervalcollection ${chromosome}_sequences_limited.fasta &
done
wait

cp all_macs_peaks_limited.narrowPeak macs_selected_chromosomes.bed
cp all_macs_peaks_limited.narrowPeak macs_selected_chromosomes_limited.bed

# Extract macs2 sequences, write to fasta (for later comparison)
echo "Extracting macs sequences"
#> macs_selected_chromosomes.bed  # Create empty file
#> macs_selected_chromosomes_limited.bed  # Create empty file
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    echo "Getting macs peaks for chromosome $chromosome"
    # Add to single file for all chromosomems
    #grep "^${chromosome}\s" $peaks_file_name >> macs_selected_chromosomes.bed
    #grep "chr${chromosome}\s" $peaks_file_name >> macs_selected_chromosomes.bed

    # Also sort out specific chromosome, making later analysis faster
    #grep "^${chromosome}\s" $peaks_file_name > macs_peaks_chr${chromosome}.bed
    #grep "chr${chromosome}\s" $peaks_file_name >> macs_peaks_chr${chromosome}.bed

    grep "^${chromosome}\s" macs_selected_chromosomes_limited.bed > macs_peaks_chr${chromosome}.bed
    grep "^chr${chromosome}\s" macs_selected_chromosomes_limited.bed >> macs_peaks_chr${chromosome}.bed

    # Compare number of macs peaks in chromosome with graph peaks. Limit to min
    #n_macs=$(wc -l macs_peaks_chr${chromosome}.bed | cut -d " " -f 1)
    #n_graph=$(wc -l ${chromosome}_max_paths.intervalcollection | cut -d " " -f 1)
    
    #echo "Chromosome $chromosome. N macs: $n_macs, N graph: $n_graph"
    #if (( n_macs > n_graph )); then
      #limit=$n_graph
    #else
      #limit=$n_macs
    #fi
    
    #fasta_limit=$(($limit*2))
    #echo "Limiting to $limit (fasta limit $fasta_limit)"

    #head -n $limit macs_peaks_chr$chromosome.bed > macs_peaks_chr${chromosome}_limited.bed
    #head -n $fasta_limit ${chromosome}_sequences.fasta > ${chromosome}_sequences_limited.fasta

    # Add macs peaks to selected chromosomes_limited
    #cat macs_peaks_chr${chromosome}_limited.bed >> macs_selected_chromosomes_limited.bed

done


# Get sequences for macs peaks
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    graph_peak_caller linear_peaks_to_fasta_summits macs_peaks_chr${chromosome}.bed $linear_genome_fasta_file  macs_sequences_chr${chromosome}.fasta $summit_window_size &

done

wait

# Fetch macs sequences for these peaks
#echo "Fetch macs sequences for selected chromosomes"
graph_peak_caller linear_peaks_to_fasta_summits macs_selected_chromosomes_limited.bed $linear_genome_fasta_file  macs_sequences.fasta $summit_window_size


# Merge all graph peak caller result files into one single sorted sequence file
graph_peak_caller concatenate_sequence_files -f "[chrom]_sequences_limited.fasta" $chromosomes sequence_all_chromosomes.fasta


# Find summits of peaks
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    graph_peak_caller get_summits -g $data_dir/$chromosome.nobg ${chromosome}_sequences_limited.fasta ${chromosome}_fragment_pileup $summit_window_size &
done

wait

graph_peak_caller concatenate_sequence_files -f "[chrom]_sequences_limited_summits.fasta" $chromosomes sequence_all_chromosomes_summits.fasta

head -n 100 sequence_all_chromosomes_summits.intervalcollection > sequence_all_chromosomes_summits_best.intervalcollection


# Project graph peaks to linear
#for chromosome in $(echo "1,2,3,4,5" | tr "," "\n")
#do
#    grep ">" ${chromosome}_sequences_summits.fasta | cut -d " " -f2-100000 > ${chromosome}_summits.intervalcollection
#	graph_peak_caller peaks_to_linear ${chromosome}_summits.intervalcollection /data/bioinf/tair2/${chromosome}_linear_pathv2.interval $chromosome ${chromosome}_peaks_to_linear.bed
#done

# Fetch sequences for projected peaks
#cat *_peaks_to_linear.bed | sort -r -n -k 9 > peaks_to_linear.bed
#graph_peak_caller linear_peaks_to_fasta peaks_to_linear.bed /data/bioinf/tair2/reference1-5.fa peaks_to_linear_sequences.fasta

#fimo --bgfile $data_dir/background.model -oc fimo_peaks_to_linear_sequences motif.meme peaks_to_linear_sequences.fasta
#$base_dir/plot_motif_enrichments.sh peaks_to_linear_sequences.fasta macs_sequences_summits.fasta $motif_url motif_enrichment_peaks_to_linear.png $tf $background_model_file


# Run motif enrichment analysis
$base_dir/plot_motif_enrichments.sh sequence_all_chromosomes_summits.fasta macs_sequences_summits.fasta $motif_url motif_enrichment.png $tf $background_model_file
cp motif_enrichment.png ../../../figures_tables/$tf.png

# Also run fimo for each chromosome
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    echo ""
    echo "----- Running fimo separately for chr $chromosome --- "
    fimo --bgfile $background_model_file -oc fimo_macs_chr$chromosome motif.meme macs_sequences_chr${chromosome}_summits.fasta
    fimo --bgfile $background_model_file -oc fimo_graph_chr$chromosome motif.meme ${chromosome}_sequences_limited_summits.fasta
done

# Analyse peak results
graph_peak_caller analyse_peaks_whole_genome --use_graph_fasta "[chrom]_sequences_limited_summits.fasta" $chromosomes ./ $data_dir ../../../figures_tables/$tf


# Hack to also run analysis for unique peaks
echo "Running analysis for unique peaks"
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    echo "Running peaks to fasta for tf $tf"
    graph_peak_caller peaks_to_fasta $data_dir/$chromosome.nobg.sequences macs_sequences_chr${chromosome}_summits_unique.intervalcollection  macs_sequences_chr${chromosome}_summits_unique.fasta
    echo "Running peaks to fasta for tf $tf (graph)"
    graph_peak_caller peaks_to_fasta $data_dir/$chromosome.nobg.sequences ${chromosome}_sequences_limited_summits_unique.intervalcollection ${chromosome}_sequences_summits_unique.fasta

    fimo --bgfile $background_model_file -oc fimo_macs_unique_chr$chromosome motif.meme macs_sequences_chr${chromosome}_summits_unique.fasta
    fimo --bgfile $background_model_file -oc fimo_graph_unique_chr$chromosome motif.meme ${chromosome}_sequences_summits_unique.fasta
done

graph_peak_caller concatenate_sequence_files -f macs_sequences_chr[chrom]_summits_unique.fasta $chromosomes unique_macs.fasta
graph_peak_caller concatenate_sequence_files -f [chrom]_sequences_summits_unique.fasta $chromosomes unique_graph.fasta

head -n 1200 unique_macs.fasta > unique_macs_top600.fasta
head -n 1200 unique_graph.fasta > unique_graph_top600.fasta

fimo --bgfile $background_model_file -oc unique_graph motif.meme unique_graph.fasta
fimo --bgfile $background_model_file -oc unique_macs motif.meme unique_macs.fasta

$base_dir/plot_motif_enrichments.sh unique_graph.fasta unique_macs.fasta $motif_url motif_enrichment_unique_peaks.png $tf $background_model_file
cp motif_enrichment_unique_peaks.png ../../../figures_tables/${tf}_unique_peaks.png

# Haplotype analysis
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    echo ""
    echo "------- Checking haplotypes --------"
    #graph_peak_caller motif_locations $data_dir ./ $chromosome
    #graph_peak_caller check_haplotype $data_dir ${data_dir}/reference1-5.fa ./ $chromosome motif_paths
done
    
