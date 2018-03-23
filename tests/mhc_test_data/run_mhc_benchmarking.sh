#!/usr/bin/env bash

graph_peak_caller create_ob_graph -o graph.nobg haplo1kg50-mhc.json

# Call peaks
graph_peak_caller callpeaks -g graph.nobg --sample reads.json -n test_ -f 135 -r 36

# Run fimo
fimo -oc fimo_test_sequences MA0139.1.meme test_sequences.fasta
fimo -oc fimo_macs_sequences MA0139.1.meme macs_sequences_mhc.fasta

# Motif enrichment
graph_peak_caller plot_motif_enrichment test_sequences.fasta macs_sequences_mhc.fasta MA0139.1.meme plot.png True CTCF

# graph_peak_caller analyse_peaks -g graph.nobg haplo1kg50-mhc.json macs_sequences_mhc.fasta test_sequences.fasta fimo_test_sequences/fimo.txt fimo_macs_sequences_mhc/fimo.txt chr6 28510119 33480577


#ls plot.png
