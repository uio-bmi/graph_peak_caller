from offsetbasedgraph import GraphWithReversals, Block, DirectedInterval
from graph_peak_caller.extender import Extender
from graph_peak_caller.sparsepileup import SparsePileup
from graph_peak_caller.areas import ValuedAreas
import logging
import numpy as np
graph = GraphWithReversals({
            i: Block(30) for i in range(1, 10010)
        },
        {i: [i+1] for i in range(1, 10009)}
    )

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

    for area in areas_list:
        valued_areas.add_binary_areas(area)

    pileup = SparsePileup.from_valued_areas(
        graph, valued_areas)

    #print(pileup)

if __name__ == "__main__":

    intervals = get_random_intervals(graph)
    for i in range(50):
        run_create_sample_pileup(intervals)