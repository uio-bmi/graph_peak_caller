#!/usr/bin/env bash

# This script runs on jenkins

./simple_chip_seq_pipeline.sh ENCSR471GSA 1 DM_JRA chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/ &
./simple_chip_seq_pipeline.sh ENCSR923VWW 1 DM_SQZ chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/ > dm_sqz_output.txt 2 >&1 &
./simple_chip_seq_pipeline.sh ENCSR978WED 1 DM_JIM chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/ &
./simple_chip_seq_pipeline.sh ENCSR082RBU 1 DM_ANTP chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/ > dm_antp_output.txt 2 >&1 &

wait

./analyse_peak_calling_results.sh  ENCSR471GSA 1 DM_JRA chr3R,chr3L,chr2R,chr2L,chrX,chr4 https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_JRA.meme /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &
./analyse_peak_calling_results.sh  ENCSR923VWW 1 DM_SQZ chr3R,chr3L,chr2R,chr2L,chrX,chr4  https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_SQZ.meme /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &
./analyse_peak_calling_results.sh  ENCSR082RBU 1 DM_ANTP chr3R,chr3L,chr2R,chr2L,chrX,chr4  https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_ANTP.meme /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &
./analyse_peak_calling_results.sh  ENCSR978WED 1 DM_JIM chr3R,chr3L,chr2R,chr2L,chrX,chr4  https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/DM_JIM.meme /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &

wait

wget -O figures_tables/bootstrap.min.css https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css
python3 generate_jenkins_report.py DM_JRA,DM_SQZ,DM_JIM,DM_ANTP figures_tables/index.html
