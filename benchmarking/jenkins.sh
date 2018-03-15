#!/usr/bin/env bash

arabidopsis_graph_dir=$1

# This script runs on jenkins

# Arabidopsis (available at NCBI: svp, ataf1, erf115, WRKY33)
# Broken fastq file. ./simple_chip_seq_pipeline.sh SRR354190 1 ARABIDOPSIS_SVP 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
# Fastq has no quality./simple_chip_seq_pipeline.sh SRX218796 1 ARABIDOPSIS_ATAF1 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/
./simple_chip_seq_pipeline.sh SRR931836 1 ARABIDOPSIS_ERF115 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRR1042995 1 ARABIDOPSIS_SEP3 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRX387187 1 ARABIDOPSIS_AP1 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRR822350	 1 ARABIDOPSIS_SOC1 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
#./simple_chip_seq_pipeline.sh SRX159033	 1 ARABIDOPSIS_PIF3 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRX310193	 1 ARABIDOPSIS_KAN1 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRR1044950 1 ARABIDOPSIS_HBI1 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &
./simple_chip_seq_pipeline.sh SRR502859 1 ARABIDOPSIS_PI 1,2,3,4,5 135000000 $arabidopsis_graph_dir/reference1-5.fa $arabidopsis_graph_dir/wg.xg $arabidopsis_graph_dir/wg.gcsa $arabidopsis_graph_dir/ &


wait

./analyse_peak_calling_results.sh  SRR931836 1 ARABIDOPSIS_ERF115 1,2,3,4,5 https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/ARABIDOPSIS_ERF115.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
./analyse_peak_calling_results.sh  SRR1042995 1 ARABIDOPSIS_SEP3 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA0563.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
./analyse_peak_calling_results.sh  SRX387187 1 ARABIDOPSIS_AP1 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA0940.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
./analyse_peak_calling_results.sh  SRR822350 1 ARABIDOPSIS_SOC1 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA0554.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
# Add later: ./analyse_peak_calling_results.sh  SRX310193 1 ARABIDOPSIS_KAN1 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA1027.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
./analyse_peak_calling_results.sh  SRR502859 1 ARABIDOPSIS_PI 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA0559.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &

# Not motif enriched: ./analyse_peak_calling_results.sh  SRX159033 1 ARABIDOPSIS_PIF3 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA0560.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
#maybe buggy? ./analyse_peak_calling_results.sh  SRR1044950 1 ARABIDOPSIS_HBI1 1,2,3,4,5 http://jaspar.genereg.net/api/v1/matrix/MA1025.1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &
#fastq problem ./analyse_peak_calling_results.sh  SRX218796 1 ARABIDOPSIS_ATAF1 1,2,3,4,5 https://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/motifs/ARABIDOPSIS_ATAF1.meme $arabidopsis_graph_dir/ $arabidopsis_graph_dir/reference1-5.fa 135000000 &

wait

wget -O figures_tables/bootstrap.min.css https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css
python3 generate_jenkins_report.py ARABIDOPSIS_ERF115,ARABIDOPSIS_SEP3,ARABIDOPSIS_AP1,ARABIDOPSIS_SOC1,ARABIDOPSIS_PI figures_tables/index.html
