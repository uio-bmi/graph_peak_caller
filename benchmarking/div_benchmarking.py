from offsetbasedgraph import GraphWithReversals, Block, DirectedInterval
from graph_peak_caller.extender import Extender
from graph_peak_caller.sparsepileup import SparsePileup, ValuedIndexes, SparseControlSample
from graph_peak_caller.areas import ValuedAreas
from graph_peak_caller.densepileup import DensePileup
import logging
import sys
import numpy as np
from graph_peak_caller.pileupcleaner2 import PeaksCleaner
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s, %(levelname)s: %(message)s")

graph = GraphWithReversals({
            i: Block(30) for i in range(1, 10010)
        },
        {i: [i+1] for i in range(1, 10009)}
    )

single_block_graph = GraphWithReversals({1: Block(30)}, {})

pileup1 = SparsePileup(single_block_graph)
pileup1.data[1] =  ValuedIndexes([4, 7, 15], [2, 5, 3], 1, 30)
pileup2 = SparsePileup(single_block_graph)
pileup2.data[1] =  ValuedIndexes([1, 10, 17], [1, 3, 7], 2, 30)

def get_random_areas(graph):
    areas = {}
    for i in range(0, 200):
        rp = np.random.randint(1, 2000)

        areas[rp] = [np.random.randint(0, 10), np.random.randint(25, 31)]

    return areas

def get_random_intervals(graph):
    intervals = []
    for i in range(0, 1500):
        rp = np.random.randint(1, 2000)

        intervals.append(
            DirectedInterval(np.random.randint(0, 10), np.random.randint(1, 10),
                     [rp, rp+1], graph)
        )
    return intervals


def run_create_sample_pileup(intervals):
    extender = Extender(graph, 80)
    valued_areas = ValuedAreas(graph)


    areas_list = (extender.extend_interval(interval)
                  for interval in intervals)
    touched_nodes = set()  # Speedup thing, keep track of nodes where areas are on
    for area in areas_list:
        valued_areas.add_binary_areas(area, touched_nodes)

    pileup = SparsePileup.from_valued_areas(
        graph, valued_areas, touched_nodes)

    #print(pileup)


def run_create_sample_pileup_new(intervals):
    extender = Extender(graph, 80)
    valued_areas = ValuedAreas(graph)
    logging.info("Extending sample reads")
    areas_list = (extender.extend_interval(interval)
                  for interval in intervals)

    logging.info("Processing areas")

    pileup = DensePileup.create_from_binary_continous_areas(
                graph, areas_list)


def graph_index():
    from graph_peak_caller.densepileupindex import GraphIndex
    graph = GraphWithReversals.from_numpy_files("../tests/lrc_kir/graph")
    index = GraphIndex.create_from_graph(graph, 200)
    index.to_file("lrc_kir")

def combine_sparsepileups(pileup1, pileup2):
    pileup = SparseControlSample.from_sparse_control_and_sample(pileup1, pileup2)


def run_peakscleaner(graph, pileup):
    cleaner = PeaksCleaner(pileup, 30)



if __name__ == "__main__":

    graph_index()
    """
    intervals = get_random_intervals(graph)
    for i in range(1):
        run_create_sample_pileup_new(intervals)
    """

    """
    pileup = SparsePileup.from_intervals(graph, get_random_intervals(graph))
    for i in range(2000):
        run_peakscleaner(graph, pileup)
    """
    """
    for i in range(10000):
        combine_sparsepileups(pileup1, pileup2)
    """

