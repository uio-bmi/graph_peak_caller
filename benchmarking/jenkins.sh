#!/usr/bin/env bash

# This script runs on jenkins

./simple_chip_seq_pipeline.sh ENCSR471GSA 1 DM_JRA chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/
#./simple_chip_seq_pipeline.sh ENCSR923VWW 1 DM_SQZ chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/ > dm_sqz_output.txt 2 >&1 &
#./simple_chip_seq_pipeline.sh ENCSR978WED 1 DM_JIM chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/ > dm_jim_output.txt 2 >&1 &
#./simple_chip_seq_pipeline.sh ENCSR082RBU 1 DM_ANTP chr3R,chr3L,chr2R,chr2L,chrX,chr4 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/ > dm_antp_output.txt 2 >&1 &

#wait

./analyse_peak_calling_results.sh  ENCSR471GSA 1 DM_JRA chr3R,chr3L,chr2R,chr2L,chrX,chr4 . /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &
#./analyse_peak_calling_results.sh  ENCSR923VWW 1 DM_SQZ chr3R,chr3L,chr2R,chr2L,chrX,chr4 . /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &
#./analyse_peak_calling_results.sh  ENCSR082RBU 1 DM_ANTP chr3R,chr3L,chr2R,chr2L,chrX,chr4 . /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &
#./analyse_peak_calling_results.sh  ENCSR978WED 1 DM_JIM chr3R,chr3L,chr2R,chr2L,chrX,chr4 . /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &

#wait
