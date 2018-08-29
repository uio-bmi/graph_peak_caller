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





# Extract macs2 sequences, write to fasta (for later comparison)
echo "Extracting macs sequences"
> macs_selected_chromosomes.bed  # Create empty file
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    echo "Getting macs peaks for chromosome $chromosome"
    # Add to single file for all chromosomems
    grep "^${chromosome}\s" $peaks_file_name >> macs_selected_chromosomes.bed
    grep "chr${chromosome}\s" $peaks_file_name >> macs_selected_chromosomes.bed

    # Also sort out specific chromosome, making later analysis faster
    grep "^${chromosome}\s" $peaks_file_name > macs_peaks_chr${chromosome}.bed
    grep "chr${chromosome}\s" $peaks_file_name >> macs_peaks_chr${chromosome}.bed
    graph_peak_caller linear_peaks_to_fasta_summits macs_peaks_chr${chromosome}.bed $linear_genome_fasta_file  macs_sequences_chr${chromosome}.fasta $summit_window_size

done

# Fetch macs sequences for these peaks
echo "Fetch macs sequences for selected chromosomes"
graph_peak_caller linear_peaks_to_fasta_summits macs_selected_chromosomes.bed $linear_genome_fasta_file  macs_sequences.fasta $summit_window_size


# Merge all graph peak caller result files into one single sorted sequence file
graph_peak_caller concatenate_sequence_files $chromosomes sequence_all_chromosomes.fasta

# Find summits of peaks
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    #graph_peak_caller get_summits -g $data_dir/$chromosome.nobg ${chromosome}_sequences.fasta ${chromosome}_qvalues $summit_window_size
    graph_peak_caller get_summits -g $data_dir/$chromosome.nobg ${chromosome}_sequences.fasta ${chromosome}_fragment_pileup $summit_window_size
done

graph_peak_caller concatenate_sequence_files -s True $chromosomes sequence_all_chromosomes_summits.fasta

head -n 100 sequence_all_chromosomes_summits.intervalcollection > sequence_all_chromosomes_summits_best.intervalcollection


# Project graph peaks to linear
for chromosome in $(echo "1,2,3,4,5" | tr "," "\n")
do
    grep ">" ${chromosome}_sequences_summits.fasta | cut -d " " -f2-100000 > ${chromosome}_summits.intervalcollection
	graph_peak_caller peaks_to_linear ${chromosome}_summits.intervalcollection /data/bioinf/tair2/${chromosome}_linear_pathv2.interval $chromosome ${chromosome}_peaks_to_linear.bed
done

# Fetch sequences for projected peaks
cat *_peaks_to_linear.bed | sort -r -n -k 9 > peaks_to_linear.bed
graph_peak_caller linear_peaks_to_fasta peaks_to_linear.bed /data/bioinf/tair2/reference1-5.fa peaks_to_linear_sequences.fasta

fimo -oc fimo_peaks_to_linear_sequences motif.meme peaks_to_linear_sequences.fasta
$base_dir/plot_motif_enrichments.sh peaks_to_linear_sequences.fasta macs_sequences_summits.fasta $motif_url motif_enrichment_peaks_to_linear.png $tf


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


# Hack to also run analysis for unique peaks
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    graph_peak_caller peaks_to_fasta $data_dir/$chromosome.nobg.sequences macs_sequences_chr${chromosome}_summits_unique.intervalcollection  macs_sequences_chr${chromosome}_summits_unique.fasta
    graph_peak_caller peaks_to_fasta $data_dir/$chromosome.nobg.sequences ${chromosome}_sequences_summits_unique.intervalcollection ${chromosome}_sequences_summits_unique.fasta

    fimo -oc fimo_macs_unique_chr$chromosome motif.meme macs_sequences_chr${chromosome}_summits_unique.fasta
    fimo -oc fimo_graph_unique_chr$chromosome motif.meme ${chromosome}_sequences_summits_unique.fasta
done

graph_peak_caller concatenate_sequence_files -f macs_sequences_chr[chrom]_summits_unique.fasta $chromosomes unique_macs.fasta
graph_peak_caller concatenate_sequence_files -f [chrom]_sequences_summits_unique.fasta $chromosomes unique_graph.fasta

head -n 1200 unique_macs.fasta > unique_macs_top600.fasta
head -n 1200 unique_graph.fasta > unique_graph_top600.fasta

fimo -oc unique_graph motif.meme unique_graph.fasta
fimo -oc unique_macs motif.meme unique_macs.fasta

$base_dir/plot_motif_enrichments.sh unique_graph.fasta unique_macs.fasta $motif_url motif_enrichment_unique_peaks.png $tf
cp motif_enrichment_unique_peaks.png ../../../figures_tables/${tf}_unique_peaks.png

# Haplotype analysis
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    echo ""
    echo "------- Checking haplotypes --------"
    #graph_peak_caller motif_locations $data_dir ./ $chromosome
    #graph_peak_caller check_haplotype $data_dir ${data_dir}/reference1-5.fa ./ $chromosome motif_paths
done
    
