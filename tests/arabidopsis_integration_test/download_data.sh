wget http://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/arabidopsis_sample.tar.gz -O sample_data.tar.gz
tar -xvzf sample_data.tar.gz

wget http://hyperbrowser.uio.no/graph-peak-caller/static/graph_peak_caller_data/tairv2.tar.gz -O graphs.tar.gz
tar -xvzf graphs.tar.gz
mv tairv2 graphs

wget http://jaspar.genereg.net/api/v1/matrix/MA0563.1.meme -O motif.meme --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6"