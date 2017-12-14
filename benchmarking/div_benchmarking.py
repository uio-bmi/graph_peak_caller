from offsetbasedgraph import GraphWithReversals, Block, DirectedInterval
from graph_peak_caller.extender import Extender
from graph_peak_caller.sparsepileup import SparsePileup, ValuedIndexes, SparseControlSample
from graph_peak_caller.areas import ValuedAreas
import logging
import numpy as np
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


def get_random_intervals(graph):
    intervals = []
    for i in range(0, 200):
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

def combine_sparsepileups(pileup1, pileup2):
    pileup = SparseControlSample.from_sparse_control_and_sample(pileup1, pileup2)



if __name__ == "__main__":

    """
    intervals = get_random_intervals(graph)
    for i in range(50):
        run_create_sample_pileup(intervals)
    """
    for i in range(10000):
        combine_sparsepileups(pileup1, pileup2)


