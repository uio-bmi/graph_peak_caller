#!/usr/bin/env bash

# This is a simple chip seq pipeline used to benchmark Graph Peak Caller.
# It run a Graph-based pipeline (using vg + graph peak caller) and a linear
# pipeline (using bwa aln + macs2) on an ENCODE Experiment.

#if (( $# != 9 )); then
#  echo "Required input arguments are missing. Check the documentation found in the Github repository."
#  exit 1
#fi

set -e

experiment_id=$1
replicate_number=$2
tf_name=$3
chromosomes=$4
genome_size=$5
linear_genome_fasta_file=$6
vg_xg_index=$7
vg_gcsa_index=$8
graph_dir=$9

base_dir=$(pwd)
n_threads=$(grep -c ^processor /proc/cpuinfo)

echo "Will fetch data from ENCODE experiment $experiment_id using replicate $replicate_number"
echo "Will use linear genome: $linear_genome_fasta_file"
echo "VG XG index: $vg_xg_index"
echo "vg gcsa index: $vg_gcsa_index"
echo "Will fetch graphs from $graph_dir"

work_dir="data/${tf_name}_${experiment_id}/$replicate_number"
mkdir -p $work_dir
cd $work_dir
echo "Made working directory $work_dir"

# Step 1: Download data
if [ ! -f raw_trimmed.fq ]; then
    echo "Will downlaod fastq"
    if [[ $experiment_id = *"ENC"* ]]; then
        echo "Experiment id is from ENCODE. Will download from ENCODE."
        encode_url=$(python3 $base_dir/download_encode_fastq.py $experiment_id $replicate_number)
        echo "Encode url: $encode_url"
        wget -O raw.fastq.gz -nv $encode_url
        echo "Unzipping"
        gunzip -c raw.fastq.gz > raw.fastq
    else
        echo "Experiment id is not from encode. Will try NCBI."
        fastq-dump $experiment_id
        mv $experiment_id.fastq raw.fastq
    fi
else
    echo "Raw fastq already exists. Not dowloading"
fi

# Step 2: Filter reads
# fastqc, trim_galore
if [ ! -f raw_trimmed.fq ]; then
    echo "Trimming fastq using trim galore"
   #mv raw.fastq filtered.fastq
   trim_galore raw.fastq
   rm raw.fastq
else
    echo "Fastq already trimmed"
fi

# Map linear reads (in order to remove bad linear mappings)
# Prepare linear reads for linear peak calling
if [ ! -f linear_alignments.bam ]; then

    echo "Mapping reads to linear genome"
    bwa aln -t $n_threads $linear_genome_fasta_file raw_trimmed.fq > reads.sai
    bwa samse $linear_genome_fasta_file  reads.sai raw_trimmed.fq > alignments.sam

    # Convert to bam and sort
    echo "Converting to bam and filtering"
    samtools view -Su alignments.sam | samtools sort - alignments_sorted

    # Filter (removed duplicates and reads having low score)
    samtools view -F 1804 -q 37 -b alignments_sorted.bam > linear_alignments.bam
fi

# Run macs with encode linear mapped reads
if [ ! -f macs_output_whole_run.txt ]; then
        echo "Running macs2"
        macs2 callpeak -m 3 100 -g $genome_size -t linear_alignments.bam -n macs > macs_output_whole_run.txt 2>&1
else
        echo "Not running max. Already done"
fi


read_length=$(cat macs_output_whole_run.txt | gawk 'match($0,  /tag size = ([0-9]+)/, ary) {print ary[1]}' )
echo "Found read length: $read_length"
fragment_length=$(cat macs_output_whole_run.txt | gawk 'match($0,  /fragment length is ([0-9]+)/, ary) {print ary[1]}' )

#echo "Found fragment length: $fragment_length"


# Step 3: Map reads
echo "Mapping reads"
if [ ! -f mapped.gam ]; then
    echo "Using indices: $vg_gcsa_index and $vg_xg_index"
    vg map -c 5000 -f raw_trimmed.fq -g $vg_gcsa_index -x $vg_xg_index > mapped.gam
    #vg mpmap --mq-method 2 -S -x $vg_xg_index -g $vg_gcsa_index -f raw_trimmed.fq > mapped.gam
else
    echo "Mapped reads exist. Not mapping"
fi

# Step 4: Filter mapped reads
echo "Filtering"
if [ ! -f filtered_low_qual_reads_removed.json ]; then
	vg filter -q 37 -t 20 mapped.gam > filtered.gam
	vg view -aj filtered.gam > filtered_low_qual_reads_removed.json
else
	echo "Filtered exists. Not filtering"
fi


# Get sequence id of reads that mapped with low mapq to linear genome
#echo "Finding reads that mapped bad to linear"
#awk '$5 < 37 { print $1  }' alignments.sam > low_qual.txt
#awk '$5 < 37 && $6 == "36M" { print $1  }' alignments.sam > low_qual.txt
#if [ ! -s "low_qual.txt" ]
#then
   #echo "Something is probaly wrong. Found no low qual reads."
   #exit 0
#fi

echo "Removing low quality reads."
#python3 $base_dir/filter_json_alignments.py low_qual.txt filtered.json filtered_low_qual_reads_removed.json
#cp filtered.json filtered_low_qual_reads_removed.json


# Get fragment length

#if [ ! -f shift_estimation_log.txt ]; then
#    graph_peak_caller estimate_shift $chromosomes $graph_dir/ filtered_low_qual_reads_removed_ 3 100 > shift_estimation_log.txt 2>&1
#else
#    echo "Not finding fragment length. Already done."
#fi
#fragment_length=$(cat shift_estimation_log.txt | gawk 'match($0,  /Found shift: ([0-9]+)/, ary) {print ary[1]}' )
#echo "Fragment length is $fragment_length"


# Step 5: Split filtered into chromosomes
if ls filtered_low_qual_reads_removed_* 1> /dev/null 2>&1; then
	echo "Not splitting into chromosomes, done before."
else
    echo "Splitting reads by chromosome."
	graph_peak_caller split_vg_json_reads_into_chromosomes $chromosomes filtered_low_qual_reads_removed.json $graph_dir
fi

# Project filtered reads onto reference
for chromosome in $(echo $chromosomes | tr "," "\n")
do 
    echo "Projecting alignments for chrom $chromosome"
    graph_peak_caller project_vg_alignments -g $graph_dir/$chromosome.nobg filtered_low_qual_reads_removed_$chromosome.json $graph_dir/${chromosome}_linear_pathv2.interval $chromosome projected_alignments_$chromosome.bed & 
done

wait  # Wait for all projections
cat projected_alignments_*.bed >> projected_alignments.bed

# Run MACS2 with projected alignments
echo "Running macs2 with projected alignments."
macs2 callpeak --nomodel --extsize $fragment_length -g $genome_size -t projected_alignments.bed -n macs_projected_alignments > macs_output_projected_alignments.txt 2>&1


# Count unique reads in filtered files
if [ ! -f count_unique_reads_output.txt ]; then
    graph_peak_caller count_unique_reads $chromosomes $graph_dir/ filtered_low_qual_reads_removed_ > count_unique_reads_output.txt 2>&1
else
    echo "Unique reads already counted. Not counting"
fi

unique_reads=$(tail -n 1 count_unique_reads_output.txt)


# Step 6 run peak caller to get p-values for every chromosome
pids=""
RESULT=0
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    if [ ! -f ${chromosome}_pvalues_values.npy ]; then
        graph_peak_caller callpeaks \
            -g $graph_dir/$chromosome.nobg \
            -s filtered_low_qual_reads_removed_$chromosome.json \
            -n ${chromosome}_ \
            -f $fragment_length \
            -r $read_length \
            -p True \
            -u $unique_reads \
            -G $genome_size \
            > log_before_p_values_$chromosome.txt 2>&1 &
            pids="$pids $!"
        echo "Peak calling (until p-values) for chr $chromosome started as process. Log will be written to $work_dir/log_before_p_values_$chromosome.txt"
    else
        echo "P values already computed for chromosome $chromosome."
    fi
done

# Wait for all to finish between continuing
for pid in $pids; do
    wait $pid || let "RESULT=1"
done

if [ "$RESULT" == "1" ];
    then
       echo "Fail in one of chromosomes. Stopping."
       exit 1
fi


# Step 7 run from p values
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    if [ -f ${chromosome}_max_paths.intervalcollection ]; then
    	echo "Peaks already called for $chromosome. Not calling"
    elif [ -f ${chromosome}_pvalues_values.npy ]; then
        graph_peak_caller callpeaks_whole_genome_from_p_values $chromosome \
            -d $graph_dir -f $fragment_length -r $read_length > log_after_p_values_$chromosome.txt 2>&1 &
	echo "Peak calling from p-values for chr $chromosome started as process. Log will be written to $work_dir/log_after_p_values_$chromosome.txt"
    else
        echo "P values not computed for $chromosome. Will not call peaks now."
    fi
done

wait
