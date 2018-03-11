
#macs2 callpeak -g 51304566 -t linear_alignments_chr22.bam -n macs -m 5 50 

graph_peak_caller callpeaks -g 22.nobg -m linear_map_22 -s reads_moved_to_graph_22.intervalcollection -n 22_ -f 146 -r 35


fimo -oc fimo_macs_chr22 motif.meme macs_sequences_chr22.fasta
fimo -oc fimo_graph_chr22 motif.meme 22_sequences.fasta

# Analyse peaks
graph_peak_caller analyse_peaks_whole_genome 22 ./ ./ results

