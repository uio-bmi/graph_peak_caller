#!/usr/bin/env bash

dm_graph_dir="/media/storage1/dm6/"
human_graph_dir="/home/ivar/data/whole_genome/"

# This script runs on jenkins

./simple_chip_seq_pipeline.sh ENCSR471GSA 1 DM_JRA chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 $dm_graph_dir/dm3_main_chromosomes.fasta $dm_graph_dir/wg.xg $dm_graph_dir/wg.gcsa $dm_graph_dir/ &
./simple_chip_seq_pipeline.sh ENCSR923VWW 1 DM_SQZ chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 $dm_graph_dir/dm3_main_chromosomes.fasta $dm_graph_dir/wg.xg $dm_graph_dir/wg.gcsa $dm_graph_dir/ > dm_sqz_output.txt 2 >&1 &
./simple_chip_seq_pipeline.sh ENCSR978WED 1 DM_JIM chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 $dm_graph_dir/dm3_main_chromosomes.fasta $dm_graph_dir/wg.xg $dm_graph_dir/wg.gcsa $dm_graph_dir/ &
./simple_chip_seq_pipeline.sh ENCSR082RBU 1 DM_ANTP chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 $dm_graph_dir/dm3_main_chromosomes.fasta $dm_graph_dir/wg.xg $dm_graph_dir/wg.gcsa $dm_graph_dir/ > dm_antp_output.txt 2 >&1 &
#./simple_chip_seq_pipeline.sh ENCSR000DUB 1 HUMAN_CTCF 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X 3080000000 $human_graph_dir/hg19_chr1-Y.fa $human_graph_dir/wg1.6.xg $human_graph_dir/wg1.6.gcsa $human_graph_dir/ > ctcf_output.txt 2 >&1 &
#./simple_chip_seq_pipeline.sh ENCSR000BIV 1 HUMAN_SRF 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X 3080000000 $human_graph_dir/hg19_chr1-Y.fa $human_graph_dir/wg1.6.xg $human_graph_dir/wg1.6.gcsa $human_graph_dir/ > srf_output.txt 2 >&1 &

wait

./analyse_peak_calling_results.sh  ENCSR471GSA 1 DM_JRA chr3R,chr3L,chr2R,chr2L,chrX,chr4 https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_JRA.meme $dm_graph_dir/ $dm_graph_dir/dm3_main_chromosomes.fasta 97958719 &
./analyse_peak_calling_results.sh  ENCSR923VWW 1 DM_SQZ chr3R,chr3L,chr2R,chr2L,chrX,chr4  https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_SQZ.meme $dm_graph_dir/ $dm_graph_dir/dm3_main_chromosomes.fasta 97958719 &
./analyse_peak_calling_results.sh  ENCSR082RBU 1 DM_ANTP chr3R,chr3L,chr2R,chr2L,chrX,chr4  https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_ANTP.meme $dm_graph_dir/ $dm_graph_dir/dm3_main_chromosomes.fasta 97958719 &
./analyse_peak_calling_results.sh  ENCSR978WED 1 DM_JIM chr3R,chr3L,chr2R,chr2L,chrX,chr4  https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_JIM.meme $dm_graph_dir/ $dm_graph_dir/dm3_main_chromosomes.fasta 97958719 &
#./analyse_peak_calling_results.sh  ENCSR000DUB 1 HUMAN_CTCF 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X  http://jaspar.genereg.net/api/v1/matrix/MA0139.1.meme $human_graph_dir/ $human_graph_dir/hg19_chr1-Y.fa 3080000000 &
#./analyse_peak_calling_results.sh  ENCSR000BIV 1 HUMAN_SRF 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X http://jaspar.genereg.net/api/v1/matrix/MA0083.2.meme $human_graph_dir/ $human_graph_dir/hg19_chr1-Y.fa 3080000000 &

wait

wget -O figures_tables/bootstrap.min.css https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css
python3 generate_jenkins_report.py DM_JRA,DM_SQZ,DM_JIM,DM_ANTP figures_tables/index.html
