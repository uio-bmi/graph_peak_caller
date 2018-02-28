import numpy as np


class SubGraphAnalyzer:
    def __init__(self, subgraph, max_path, dists):
        self._subgraph = subgraph
        self._score = np.min(dists[:, -1])
        self._max_path = max_path
        self._dists = dists

    def get_info(self):
        return (self.has_two_bindings(), self.is_ambiguous())

    def has_two_bindings(self):
        if len(self._max_path) <= 1:
            return False
        sizes = np.min(self._subgraph, 1)
        return np.sum(sizes) != self._score

    def is_ambiguous(self):
        if len(self._max_path) <= 1:
            return False
        for node in self._max_path[:-1]:
            to_dist = self._dists[self._max_path[0], node]
            row = self._subgraph[node]
            for next_node, d in zip(row.indices, row.data):
                if next_node in self._max_path:
                    continue
                s = to_dist + d + self._dists[next_node, -1]
                assert s >= self._score
                if s == self._score:
                    return True

        return False
