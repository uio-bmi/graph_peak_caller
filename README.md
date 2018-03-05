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
Graph Peak Caller is typically used together with [vg](http://github.com/vgteam/vg) and can take as input alignments produces by vg. Graph Peak Caller uses the Python module [pyvg](https://github.com/uio-bmi/pyvg) (installed together with Graph Peak Caller) to convert output from vg to formats compatible with Python.

If you have a vg graph

## Advanced usage

