import logging
from .control.linearmap import LinearMap


def create_linear_map(ob_graph, out_file_name="linear_map.npz"):
    linear_map = LinearMap.from_graph(ob_graph)
    linear_map.to_file(out_file_name)
    logging.info("Created linear map, wrote to file %s" % out_file_name)
