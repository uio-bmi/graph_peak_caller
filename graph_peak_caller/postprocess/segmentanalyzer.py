import numpy as np
from itertools import chain


class SegmentAnalyzer:
    def __init__(self, segments, node_ids, node_indexes):
        self._segments = segments
        self._node_ids = node_ids
        self.node_indexes = node_indexes
        self.fulls = []
        self.starts = []
        self.ends = []

    def _get_spanning_nodes(self, start_ids, end_ids):
        spanning = end_ids > start_ids+1
        return list(chain.from_iterable(
            range(start+1, end) for start, end in
            zip(start_ids[spanning], end_ids[spanning])))

    def _is_starts(self, offsets, node_ids):
        return offsets == self.node_indexes[node_ids-1]

    def _is_ends(self, offsets, node_ids):
        return offsets == self.node_indexes[node_ids]

    def get_starts(self):
        a = np.transpose(np.array(self.starts))
        if not a.size:
            return np.empty((2, 0), dtype="int")
        a[1] = a[1]-self.node_indexes[a[0]-1]
        return a.astype("int")

    def get_ends(self):
        a = np.transpose(np.array(self.ends, dtype="int"))
        if not a.size:
            return np.empty((2, 0), dtype="int")
        a[1] = self.node_indexes[a[0]]-a[1]
        return a.astype("int")

    def get_fulls(self):
        a = np.array(self.fulls, dtype="int")
        if not a.size:
            return np.empty((2, 0))
        sizes = self.node_indexes[a] - self.node_indexes[a-1]
        return np.vstack((a, sizes.astype("int"))).astype("int")

    def add_starts(self, node_ids, offsets):
        self.starts.extend(zip(node_ids, offsets))

    def add_ends(self, node_ids, offsets):

        self.ends.extend(zip(node_ids, offsets))

    def add_fulls(self, node_ids):
        self.fulls.extend(node_ids)

    def handle_internal(self):
        is_starts = self._is_starts(self._internal_segments[:, 0],
                                    self._internal_ids[:, 0])
        is_ends = self._is_ends(self._internal_segments[:, 1],
                                self._internal_ids[:, 1])
        full_mask = is_starts & is_ends
        self.add_fulls(self._internal_ids[:, 0][full_mask])
        start_mask = is_starts & ~full_mask
        self.add_starts(self._internal_ids[:, 1][start_mask],
                        self._internal_segments[:, 1][start_mask])
        end_mask = is_ends & ~full_mask
        self.add_ends(self._internal_ids[:, 0][end_mask],
                      self._internal_segments[:, 0][end_mask])
        self.internals = self._internal_segments[~is_starts & ~is_ends]

    def handle_multinode(self):
        self.add_fulls(self._get_spanning_nodes(self._spanning_ids[:, 0],
                                                self._spanning_ids[:, 1]))
        is_starts = self._is_starts(self._spanning_segments[:, 0],
                                    self._spanning_ids[:, 0])
        is_ends = self._is_ends(self._spanning_segments[:, 1],
                                self._spanning_ids[:, 1])
        self.add_fulls(self._spanning_ids[:, 0][is_starts])
        self.add_fulls(self._spanning_ids[:, 1][is_ends])
        self.add_ends(self._spanning_ids[:, 0][~is_starts],
                      self._spanning_segments[:, 0][~is_starts])
        self.add_starts(self._spanning_ids[:, 1][~is_ends],
                        self._spanning_segments[:, 1][~is_ends])

    def run(self):
        multinodes = self._node_ids[:, 0] != self._node_ids[:, 1]
        self._internal_segments = self._segments[~multinodes]
        self._internal_ids = self._node_ids[~multinodes]
        self._spanning_segments = self._segments[multinodes]
        self._spanning_ids = self._node_ids[multinodes]
        self.multinodes = multinodes
        self.handle_multinode()
        self.handle_internal()


class SegmentSplitter(SegmentAnalyzer):
    def handle_internal(self):
        is_starts = self._is_starts(self._internal_segments[:, 0],
                                    self._internal_ids[:, 0])
        is_ends = self._is_ends(self._internal_segments[:, 1],
                                self._internal_ids[:, 1])
        self.full_mask = is_starts & is_ends
        self.start_mask = is_starts & ~self.full_mask
        self.end_mask = is_ends & ~self.full_mask
        self.internal_mask = ~is_starts & ~is_ends

    def handle_multinode(self):
        spanned = self._get_spanning_nodes(self._spanning_ids[:, 0],
                                           self._spanning_ids[:, 1])
        n_spanned = len(spanned)
        n_old = self._segments.shape[0]
        n_new = self._spanning_segments.shape[0]
        n_total = n_spanned + n_old + n_new
        splitted_segments = np.empty((n_total, 2), dtype="int")
        if n_old:
            splitted_segments[:n_old] = self._segments
            splitted_segments[:n_old, 1][self.multinodes] = self.node_indexes[
                self._spanning_ids[:, 0]]
        splitted_segments[n_old:n_old+n_new, 0] = self.node_indexes[
            self._spanning_ids[:, 1]-1]
        splitted_segments[n_old:n_old+n_new, 1] = self._spanning_segments[:, 1]
        spanned = np.array(spanned, dtype="int")
        if n_spanned:
            splitted_segments[-n_spanned:, 0] = self.node_indexes[spanned-1]
            splitted_segments[-n_spanned:, 1] = self.node_indexes[spanned]
        self.splitted_segments = splitted_segments
        self._internal_segments = splitted_segments
        ids = np.r_[self._node_ids[:, 0],
                    self._node_ids[:, 1][self.multinodes],
                    spanned]
        # assert np.all(ids >= self._graph.blocks.node_id_offset), (self._graph.blocks.node_id_offset,
        #                                                           ids[ids<self._graph.blocks.node_id_offset])
        self._internal_ids = np.hstack((ids[:, None], ids[:, None]))


def test_segments_cleaner():
    offsets = np.array([80, 100,
                        180, 220,
                        240, 250,
                        300, 400,
                        500, 520,
                        610, 810])

    offsets = offsets.reshape((offsets.size//2, 2))
    node_ids = np.empty_like(offsets)
    node_ids[:, 0] = offsets[:, 0] // 100
    node_ids[:, 1] = ((offsets[:, 1]-1) // 100)
    node_ids += 1
    node_indexes = np.arange(0, 1201, 100)
    segments = SegmentsAnalyzer(offsets, node_ids, node_indexes)
    segments = np.unique(np.sort(segments.ravel()))

if __name__ == "__main__":
    test_segments_cleaner()
    
