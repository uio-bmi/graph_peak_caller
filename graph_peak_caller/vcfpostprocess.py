class MaxPaths:
    def __init__(self, graph, edge_pileup, snp_pileup):
        self._graph = graph
        self._edge_pileup = edge_pileup
        self._snp_pileup = snp_pileup

    def 


class PostProcess:
    def __init__(self, graph, linear_repr, max_hole_size, min_peak_length):
        self._graph = graph
        self._linear_repr = linear_repr
        self._max_hole_size = max_hole_size
        self._min_peak_length = min_peak_length

    def fill_holes(self, peaks):
        pass

    def remove_small_peaks(self, peaks):


def fill_holes(graph, linear_repr, peaks):
    # Find peak-ends/hole-starts
    # Separate into internal and spanning
    # For internal:
    #     Filter out big ones
    # For spanning:
    # Filter out those with large distance to node-end_id
    #
