#!/usr/bin/env bash

graph_peak_caller create_ob_numpy_graph haplo1kg50-mhc.json graph
graph_peak_caller create_ob_graph haplo1kg50-mhc.json graph.obg

# Create snarls
graph_peak_caller create_linear_map graph.obg haplo1kg50-mhc.snarls linear_map

# Call peaks
graph_peak_callery callpeaks_with_numpy_graph graph haplo1kg50-mhc.vg linear_map reads.json reads.json False test_ 135 36

# Motif enrichment
#graph_peak_caller plot_motif_enrichment test_sequences.fasta macs_sequences.fasta MA0139.1.meme plot.png

#ls plot.png