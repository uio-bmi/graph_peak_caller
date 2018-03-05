[![Build Status](https://travis-ci.org/uio-bmi/graph_peak_caller.svg?branch=master)](https://travis-ci.org/uio-bmi/graph_peak_caller)
[![codecov](https://codecov.io/gh/uio-bmi/graph_peak_caller/branch/master/graph/badge.svg)](https://codecov.io/gh/uio-bmi/graph_peak_caller)

# Graph Peak Caller
Graph Peak Caller is a tool for calling ChIP-seq peaks on graph-based reference genomes. Graph Peak Caller is easiest to use together with [vg](http://github.com/vgteam/vg) and can be used both as a command-line tool and as a Python module.

## Installation
Graph Peak Caller is written in Python 3 and can be installed through Pip:
```
pip3 install graph_peak_caller
```

Validate that the installation worked by writing `graph_peak_caller version` on your command line.

# Usage guide
## Use Graph Peak Caller through Galaxy
Graph Peak Caller can be used through a Galaxy installation at [hyperbrowser.uio.no/graph-peak-caller](hyperbrowser.uio.no/graph-peak-caller). A set of pre-compiled graphs are available to use (see the welcome page at the Galaxy server for an overview). If one of these graphs suites your needs, this is the best way to use Graph Peak Caller. PS: If you would like us to include graph-based reference genomes for a new species, please contact us, e.g by opening an issue, and we will do our best. 

## Using Graph Peak Caller on the command line with vg
Graph Peak Caller made to be used together with [vg](http://github.com/vgteam/vg) and take as input alignments produces by vg. Graph Peak Caller uses the Python module [pyvg](https://github.com/uio-bmi/pyvg) (installed together with Graph Peak Caller) to convert output from vg to formats compatible with Python.

If you have a single vg graph representing your reference genome, the following explains you how to use Graph Peak Caller with that graph. If you have multiple vg graphs, see *Advanced usage* below.

### Step 1: Preparing data
Convert your vg graph to json and create an Offset Based Python Graph:
```
vg view -Jv graph.vg > graph.json
graph_peak_caller create_ob_graph graph.json graph.nobg
```

Graph Peak Caller uses a "linear map" for projecting the graph structure to a linear structure when creating the background signal. This can be created from a vg snarl graph:
```
vg snarls graph.vg > snarls.pb
graph_peak_caller create_linear_map graph.nobg snarls.pb linear_map
```

Also, convert your vg alignments to json:
```
vg view -aj alignments.gam > alignments.json
vg view -aj control_alignments.gam > control_alignments.json
```
If you don't have any vg alignments, and do not know how to produce them, check out [this vg guide](https://github.com/vgteam/vg/wiki/Basic-Operations).

### Step 2: Call peaks
Finally, we can call peaks by using the *callpeaks* command:
```
graph_peak_caller callpeaks graph.nobg graph.json linear_map alignments.json control_alignments.json True test FRAGMENT_LENGTH READ_LENGTH
```

Make sure to change *FRAGMENT_LENGTH* and *READ_LENGTH* with numbers matching your data. If you do not know the fragment length of your ChIP-seq experiment, you can use Macs' predictd command to find estimate it: `macs predict -t alignments.bam`.

If you do not have a control track, use your sample reads as control and change True to False in the above command. Chaning True to False is important as Graph Peak Caller will generate the background signal in a different way when the sample is used as control. 

## Advanced usage
If you want to do Peak Calling on a whole-genome reference, vg will typically produce one graph for each chromosome, and it is best to divide the peak calling into one process for each chromosome. Check out [this guide](https://github.com/uio-bmi/graph_peak_caller/wiki/Graph-based-ChIP-seq-tutorial) for a detailed explanation on how that can be done.

# Development
