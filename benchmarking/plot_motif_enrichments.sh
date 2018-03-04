#!/usr/bin/env bash
source config_server.sh
echo "Graph peak caller: $graph_peak_caller"
echo "Plot motif enrichments"

fasta1=$1
fasta2=$2
motif_url=$3
out_file=$4
tf=$5

if [ ! -f motif.meme ]; then
	wget -O motif.meme --show-progress $motif_url --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6"
fi

graph_peak_caller plot_motif_enrichment $fasta1 $fasta2 motif.meme $out_file True $tf

