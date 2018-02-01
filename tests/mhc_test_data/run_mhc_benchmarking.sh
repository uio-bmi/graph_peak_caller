#!/usr/bin/env bash

python3 ../../graph_peak_caller.py create_ob_numpy_graph haplo1kg50-mhc.json graph
python3 ../../graph_peak_caller.py create_ob_graph haplo1kg50-mhc.json graph.obg

# Create snarls
python3 ../../graph_peak_caller.py create_linear_map graph.obg haplo1kg50-mhc.snarls linear_map

# Call peaks
python3 ../../graph_peak_caller.py callpeaks_with_numpy_graph graph haplo1kg50-mhc.vg linear_map reads.json reads.json False test_ 135 36

# Motif enrichment
#python3 ../../graph_peak_caller.py plot_motif_enrichment test_sequences.fasta macs_sequences.fasta MA0139.1.meme plot.png

#ls plot.png