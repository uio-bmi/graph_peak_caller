#!/usr/bin/env bash

arabidopsis_graph_dir=$1
dm_graph_dir=$2
human_graph_dir=$3

echo "Using graph dirs $arabidopsis_graph_dir, $dm_graph_dir and $human_graph_dir"

# This script runs on jenkins

./simple_chip_seq_pipeline.sh SRR931836 1 ARABIDOPSIS_ERF115 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRR1042995 1 ARABIDOPSIS_SEP3 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRX387187 1 ARABIDOPSIS_AP1 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRR822350	 1 ARABIDOPSIS_SOC1 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRR1044950 1 ARABIDOPSIS_HBI1 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRR502859 1 ARABIDOPSIS_PI 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &

# Drosophila
if [ ! -z "$dm_graph_dir" ]
then
    ./simple_chip_seq_pipeline.sh ENCSR471GSA 1 DM_JRA chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 $dm_graph_dir/dm3_main_chromosomes.fasta $dm_graph_dir/wg.xg $dm_graph_dir/wg.gcsa $dm_graph_dir/ &
    ./simple_chip_seq_pipeline.sh ENCSR923VWW 1 DM_SQZ chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 $dm_graph_dir/dm3_main_chromosomes.fasta $dm_graph_dir/wg.xg $dm_graph_dir/wg.gcsa $dm_graph_dir/ > dm_sqz_output.txt 2 >&1 &
    ./simple_chip_seq_pipeline.sh ENCSR978WED 1 DM_JIM chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 $dm_graph_dir/dm3_main_chromosomes.fasta $dm_graph_dir/wg.xg $dm_graph_dir/wg.gcsa $dm_graph_dir/ &
    ./simple_chip_seq_pipeline.sh ENCSR082RBU 1 DM_ANTP chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 $dm_graph_dir/dm3_main_chromosomes.fasta $dm_graph_dir/wg.xg $dm_graph_dir/wg.gcsa $dm_graph_dir/ > dm_antp_output.txt 2 >&1 &
fi

# Human
if [ ! -z "$human_graph_dir" ]
then
    ./simple_chip_seq_pipeline.sh ENCSR000DUB 1 HUMAN_CTCF 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X 3080000000 $human_graph_dir/hg19_chr1-Y.fa $human_graph_dir/wg1.6.xg $human_graph_dir/wg1.6.gcsa $human_graph_dir/ > ctcf_output.txt 2 >&1 &
    ./simple_chip_seq_pipeline.sh ENCSR000BIV 1 HUMAN_SRF 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X 3080000000 $human_graph_dir/hg19_chr1-Y.fa $human_graph_dir/wg1.6.xg $human_graph_dir/wg1.6.gcsa $human_graph_dir/ > srf_output.txt 2 >&1 &
fi

wait

./analyse_peak_calling_results.sh  SRR931836 1 ARABIDOPSIS_ERF115 1,2,3,4,5 https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/ARABIDOPSIS_ERF115.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
./analyse_peak_calling_results.sh  SRR1042995 1 ARABIDOPSIS_SEP3 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA0563.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
./analyse_peak_calling_results.sh  SRX387187 1 ARABIDOPSIS_AP1 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA0940.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
./analyse_peak_calling_results.sh  SRR822350 1 ARABIDOPSIS_SOC1 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA0554.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
./analyse_peak_calling_results.sh  SRR502859 1 ARABIDOPSIS_PI 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA0559.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &

if [ ! -z "$dm_graph_dir" ]
then
    ./analyse_peak_calling_results.sh  ENCSR471GSA 1 DM_JRA chr3R,chr3L,chr2R,chr2L,chrX,chr4 https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_JRA.meme $dm_graph_dir/ $dm_graph_dir/dm3_main_chromosomes.fasta 97958719 &
    ./analyse_peak_calling_results.sh  ENCSR923VWW 1 DM_SQZ chr3R,chr3L,chr2R,chr2L,chrX,chr4  https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_SQZ.meme $dm_graph_dir/ $dm_graph_dir/dm3_main_chromosomes.fasta 97958719 &
    ./analyse_peak_calling_results.sh  ENCSR082RBU 1 DM_ANTP chr3R,chr3L,chr2R,chr2L,chrX,chr4  https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_ANTP.meme $dm_graph_dir/ $dm_graph_dir/dm3_main_chromosomes.fasta 97958719 &
    ./analyse_peak_calling_results.sh  ENCSR978WED 1 DM_JIM chr3R,chr3L,chr2R,chr2L,chrX,chr4  https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_JIM.meme $dm_graph_dir/ $dm_graph_dir/dm3_main_chromosomes.fasta 97958719 &
fi

if [ ! -z "$human_graph_dir" ]
then
    ./analyse_peak_calling_results.sh  ENCSR000DUB 1 HUMAN_CTCF 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X  http://jaspar.genereg.net/api/v1/matrix/MA0139.1.meme $human_graph_dir/ $human_graph_dir/hg19_chr1-Y.fa 3080000000 &
    ./analyse_peak_calling_results.sh  ENCSR000BIV 1 HUMAN_SRF 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X http://jaspar.genereg.net/api/v1/matrix/MA0083.2.meme $human_graph_dir/ $human_graph_dir/hg19_chr1-Y.fa 3080000000 &
fi


wait

wget -O figures_tables/bootstrap.min.css https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css

if [! -z "${human_graph_dir}$dm_graph_dir" ]; then
    python3 generate_jenkins_report.py ARABIDOPSIS_ERF115,ARABIDOPSIS_SEP3,ARABIDOPSIS_AP1,ARABIDOPSIS_SOC1,ARABIDOPSIS_PI,DM_JRA,DM_SQZ,DM_JIM,DM_ANTP,HUMAN_CTCF,HUMAN_SRF figures_tables/index.html
else
    # Only run arabidopsis
    python3 generate_jenkins_report.py ARABIDOPSIS_ERF115,ARABIDOPSIS_SEP3,ARABIDOPSIS_AP1,ARABIDOPSIS_SOC1,ARABIDOPSIS_PI figures_tables/index.html
fi