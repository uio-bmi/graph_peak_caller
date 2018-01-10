from pyvg.util import vg_json_file_to_interval_collection, vg_gam_file_to_interval_collection
import offsetbasedgraph as obg
from graph_peak_caller.extender import Extender
from graph_peak_caller.areas import ValuedAreas
from graph_peak_caller.sparsepileupv2 import SparsePileup
from graph_peak_caller.sparsepileup import SparsePileup as OldSparsePileup
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s, %(levelname)s: %(message)s")


ob_graph = obg.GraphWithReversals.from_numpy_files("graph__")
#intervals = vg_gam_file_to_interval_collection(None, "reads.gam", ob_graph)
intervals = vg_json_file_to_interval_collection(None, "../tests/lrc_kir/small.json", ob_graph)

i = 0

def check_get_intervals_time():
    i = 0
    for interval in intervals:
        #print(interval)
        if i % 5000 == 0:
            print("Interval %d" % i)

        i += 1



def run_extender():
    extender = Extender(ob_graph, 140)
    valued_areas = ValuedAreas(ob_graph)
    areas_list = (extender.extend_interval(interval)
                  for interval in intervals)

    i = 0
    touched_nodes = set()
    for area in areas_list:
        if i % 5000 == 0:
            print("Processing area %d" % i)
        i += 1

        valued_areas.add_binary_areas(area, touched_nodes)

    pileup = SparsePileup.from_valued_areas(
            ob_graph, valued_areas, touched_nodes=touched_nodes)



run_extender()