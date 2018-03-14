#!/usr/bin/env bash

encode_id=$1
replicate=$2
tf=$3
chromosomes=$4
motif_url=$5
data_dir=$6
linear_genome_fasta_file=$7
genome_size=$8

base_dir=$(pwd)
n_threads=$(grep -c ^processor /proc/cpuinfo)

echo "Running analyse_peak_calling_results.sh"
work_dir=data/${tf}_${encode_id}/${replicate}
cd $work_dir

echo "Changed dir to $work_dir"

# Extract macs2 sequences, write to fasta (for later comparison)
echo "Extracting macs sequences"
> macs_selected_chromosomes.bed  # Create empty file
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    echo "Getting macs peaks for chromosome $chromosome"
    # Add to single file for all chromosomems
    grep "^${chromosome}\s" macs_peaks.narrowPeak >> macs_selected_chromosomes.bed
    grep "chr${chromosome}\s" macs_peaks.narrowPeak >> macs_selected_chromosomes.bed

    # Also sort out specific chromosome, making later analysis faster
    grep "^${chromosome}\s" macs_peaks.narrowPeak > macs_peaks_chr${chromosome}.bed
    grep "chr${chromosome}\s" macs_peaks.narrowPeak >> macs_peaks_chr${chromosome}.bed
    graph_peak_caller linear_peaks_to_fasta macs_peaks_chr${chromosome}.bed $linear_genome_fasta_file  macs_sequences_chr${chromosome}.fasta

done

# Fetch macs sequences for these peaks
echo "Fetch macs sequences for selected chromosomes"
graph_peak_caller linear_peaks_to_fasta macs_selected_chromosomes.bed $linear_genome_fasta_file  macs_sequences.fasta


# Merge all graph peak caller result files into one single sorted sequence file
graph_peak_caller concatenate_sequence_files $chromosomes sequence_all_chromosomes.fasta

# Find summits of peaks
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    graph_peak_caller get_summits -g graphs/$chromosome.nobg ${chromosome}_sequences.fasta ${chromosome}_qvalues
done
graph_peak_caller concatenate_sequence_files -s True $chromosomes sequence_all_chromosomes_summits.fasta

# Run motif enrichment analysis
$base_dir/plot_motif_enrichments.sh sequence_all_chromosomes_summits.fasta macs_sequences_summits.fasta $motif_url motif_enrichment.png $tf
cp motif_enrichment.png ../../../figures_tables/$tf.png

# Also run fimo for each chromosome
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    echo ""
    echo "----- Running fimo separately for chr $chromosome --- "
    fimo -oc fimo_macs_chr$chromosome motif.meme macs_sequences_chr${chromosome}_summits.fasta
    fimo -oc fimo_graph_chr$chromosome motif.meme ${chromosome}_sequences_summits.fasta
done

# Analyse peak results
graph_peak_caller analyse_peaks_whole_genome $chromosomes ./ $data_dir ../../../figures_tables/$tf
