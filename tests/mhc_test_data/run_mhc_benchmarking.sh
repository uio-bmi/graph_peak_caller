#!/usr/bin/env bash

graph_peak_caller create_ob_graph haplo1kg50-mhc.json graph.obg

# Create snarls
graph_peak_caller create_linear_map graph.obg haplo1kg50-mhc.snarls linear_map

# Call peaks
graph_peak_caller callpeaks graph.obg haplo1kg50-mhc.vg linear_map reads.json reads.json False test_ 135 36

# Run fimo
fimo -oc fimo_test_sequences MA0139.1.meme test_sequences.fasta
fimo -oc fimo_macs_sequences MA0139.1.meme macs_sequences.fasta

# Motif enrichment
graph_peak_caller plot_motif_enrichment test_sequences.fasta macs_sequences.fasta MA0139.1.meme plot.png True

# graph_peak_caller analyse_peaks graph.obg haplo1kg50-mhc.json macs_sequences_chr6.fasta test_sequences.fasta fimo_test_sequences/fimo.txt fimo_macs_sequences/fimo.txt


#ls plot.png