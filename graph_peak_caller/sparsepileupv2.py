import logging
from itertools import chain
import numpy as np
from scipy.stats import poisson
from collections import defaultdict
from .pileup import Pileup
from .pileupcleaner2 import PeaksCleaner, HolesCleaner
from .subgraphcollection import SubgraphCollection
from .eventsorter import DiscreteEventSorter
from offsetbasedgraph import Interval, IntervalCollection
import pickle
from .sparsepileup import SparseAreasDict, starts_and_ends_to_sparse_pileup, intervals_to_start_and_ends
from memory_profiler import profile

class SimpleValuedIndexes():
    def __init__(self, indexes, values):
        self.indexes = indexes
        self.values = values

    def sum(self):
        lengths = np.diff(self.indexes)
        return np.sum(lengths*self.values)


class RpScore:
    def __init__(self, max_score, sum_score):
        self.sum_score = sum_score
        self.max_score = max_score

    def __getitem__(self, item):
        if item == 0:
            return self.max_score
        elif item == 1:
            return self.sum_score
        else:
            raise NotImplementedError()

    def sum(self):
        return self.sum_score

    def max_value(self):
        return self.max_score

    @classmethod
    def from_valued_indexes(cls, vi):
        return cls(np.max(vi.all_values()), vi.sum())


class SparsePileupData:
    def __init__(self, node_ids, lengths, ndim=1, graph=None, default_value=0):

        self._node_indexes = None  # Mapping from node id to index in indexes/values
        self._indexes = None
        self._values = None
        self._lengths = None
        self.min_node = None
        self.min_value = 0
        self.default_value = default_value
        self.graph = graph
        self.nodes = set()

        if len(node_ids) > 0:
            self._create_empty(node_ids, lengths, ndim)

    def _create_empty(self, node_ids, lengths, ndim):

        self.nodes = set(node_ids)

        sorted_nodes = sorted(node_ids)
        self.min_node = sorted_nodes[0]
        max_node = sorted_nodes[-1]
        span = max_node - self.min_node + 1

        n_elements = sum(lengths)
        self._indexes = np.zeros(n_elements, dtype=np.uint16)
        if ndim > 1:
            self._values = np.zeros((n_elements, ndim), dtype=np.float16)
        else:
            self._values = np.zeros(n_elements)
        self._node_indexes = np.zeros(span, dtype=np.uint32)
        self._lengths = np.zeros(span, dtype=np.uint16)

        logging.info("N nodes in sparse pileup: %d" % len(self.nodes))

        offset = 0
        for i, node in enumerate(node_ids):
            index = node - self.min_node
            self._node_indexes[index] = offset
            self._lengths[index] = lengths[i]
            assert lengths[i] >= 2
            assert lengths[i] <= 1000
            #print("Node %d at offset %d" % (node, offset))
            offset += lengths[i]

    def copy(self):
        new = SparsePileupData([], [], graph=self.graph,
                               default_value=self.default_value)
        new._indexes = self._indexes.copy()
        new._values = self._values.copy()
        new._node_indexes= self._node_indexes.copy()
        new._lengths = self._lengths.copy()
        new.min_node = self.min_node
        new.nodes = self.nodes.copy()

        return new

    def __eq__(self, other):
        assert isinstance(other, SparsePileupData)
        for node in self.nodes:
            if not np.all(self.indexes(node) == other.indexes(node)):
                print("Indices %s != %s" % (self.indexes(node), other.indexes(node)))
                return False

            assert isinstance(other.values(node), np.ndarray)
            print("Other values")
            print(other.values(node))

            if np.all(np.abs(self.values(node) -  other.values(node)) > 1e-5):
                print("Values %s != %s" % (self.values(node), other.values(node)))
                print()
                return False

        return True

    def __len__(self):
        return len(self.nodes)

    def set_min_value(self, value):
        self.min_value = value

    def set_values(self, node_id, values):
        assert self._values is not None
        # First index is always start value
        # This method only sets everything else than start
        assert node_id in self.nodes, "Only allowed to set for already initiated nodes"
        index = node_id - self.min_node
        start = self._node_indexes[index] + 1
        end = start + len(values)
        self._values[start:end] = values

    def set_indexes(self, node_id, indexes):
        index = node_id - self.min_node
        if len(indexes) == 0:
            self._lengths[index] = 2
            return

        assert np.all(np.diff(indexes) > 0), "Invalid indexes %s" % indexes
        # First index is always 0. Last index is always length
        # This method sets everything else than first index
        assert node_id in self.nodes, "Only allowed to set for already initiated nodes"
        start = self._node_indexes[index] + 1
        end = start + len(indexes)
        assert len(indexes) <= self._lengths[index] , "Maybe allocate space for higher np.integers in _lengths"
        #assert self._indexes[end] == 0 or self._indexes[end] > indexes[-1], "Next index is %d when setting indexes %s" % (self._indexes[end], indexes)
        self._indexes[start:end] = indexes
        self._lengths[index] = len(indexes) + 2  # Can be removed if handled in
                                                 # special cases outside (when more space is allocated than needed)
        assert self._indexes[start-1] == 0, "Start index should be 0, is %d" % self._indexes[start-1]

        assert np.all(self.indexes(node_id)[1:-1] == indexes), "%s != %s" % (self.indexes(node_id)[1:-1], indexes)

        assert np.all(np.diff(self.indexes(node_id)[:-1]) > 0), "Invalid indexes %s after setting indexes %s" % (self.indexes(node_id), indexes)

    def set_end_index(self, node_id, i):
        # Used to set "length" of node. End should always be length
        index = node_id - self.min_node
        pos =  self._node_indexes[index] + self._lengths[index] - 1

        #for idx in self.indexes(node_id)[:-1]:
        #    assert i > idx, "Trying to set end idx %d which is <= already existing index %d. All indexes: %s" % (i, idx, self.indexes(node_id))

        #assert self._indexes[end-1] == 0, "End index already has something"
        self._indexes[pos] = i
        assert self._indexes[pos] == i, "%d != %d" % (self._indexes[pos], i)

        assert self.indexes(node_id)[-1] == i, "End index %d is not %d. All indexes: %s" % \
            (i, self.indexes(node_id)[-1], self.indexes(node_id))


    def set_start_value(self, node_id, val):
        index = node_id - self.min_node
        start = self._node_indexes[index]
        self._values[start] = val

    def indexes(self, node_id):
        if node_id not in self.nodes:
            return np.array([0, self.graph.node_size(node_id)])

        index = node_id - self.min_node

        length = self._lengths[index]
        assert length >= 2
        start = self._node_indexes[index]
        end = start + length
        return self._indexes[start:end]

    def values(self, node_id):
        if node_id not in self.nodes:
            return np.array([self.default_value])

        index = node_id - self.min_node
        start = self._node_indexes[index]
        end = start + self._lengths[index] - 1
        return self._values[start:end]

    def sum(self):
        assert self.min_value == 0
        return np.sum(self._values)

    def scale(self, factor):
        self._values *= factor

    def fill_existing_hole(self, node, index, value):
        # Fills an existing hole starting at index in node
        idx = np.nonzero(self.indexes(node) == index)
        assert len(idx) == 1, "There should be only one index matching the hole"
        #assert idx == len(self.indexes) -1], "Trying to set value for end index (has no value)"
        old_values = self.values(node)
        assert old_values[idx] == 0
        old_values[idx] = value
        if idx == 0:
            self.set_start_value(value)
        else:
            self.set_values(node, old_values[1:])

    def sanitize_node_indices(self, node_id):
        indexes = self.indexes(node_id)
        values = self.values(node_id)
        #print("Sanitizing")
        #print("INdexes: %s" % indexes)
        #print("Values: %s" % values)
        diffs = np.where(np.diff(indexes) > 0)[0]
        #print("New indexes: %s" % indexes[diffs])
        self.set_indexes(node_id, indexes[diffs[1:]])
        self.set_end_index(node_id, indexes[-1])
        self.set_values(node_id, values[diffs[1:]])
        #self._lengths[node_id - self.min_node] = len(diffs)

    def sanitize_node(self, node_id):
        all_vals = self.values(node_id)
        indexes = self.indexes(node_id)
        assert len(indexes) >= 2
        #print("Old indices: %s" % indexes)
        #print("Old values: %s" % all_vals)

        diffs = np.diff(all_vals)
        changes = np.where(diffs != 0)[0]
        if len(all_vals) == 0:
            return

        new_indexes = indexes[changes+1]
        new_values = all_vals[changes+1]
        #print("New indices: %s" % new_indexes)
        #print("New values: %s" % new_values)
        self.set_values(node_id, new_values)
        self.set_indexes(node_id, new_indexes)
        self.set_end_index(node_id, indexes[-1])

        new_indexes = self.indexes(node_id)
        assert len(new_indexes) >= 2
        assert new_indexes[0] == 0
        assert new_indexes[-1] == indexes[-1]

    def get_flat_numpy_array(self, node_id):
        indexes = self.indexes(node_id)
        flat = np.zeros(indexes[-1])
        for s, e, value in self.index_value_pairs(node_id):
            flat[s:e] = value
        return flat

    def get_subset_max_value(self, node_id, start, end):
        return np.max(self.get_flat_numpy_array(node_id)[start:end])

    def get_subset_sum(self, node_id, start, end):
        return np.sum(self.get_flat_numpy_array(node_id)[start:end])

    def score(self, node_id, start, end):
        return RpScore(self.get_subset_max_value(node_id, start, end), \
               self.get_subset_sum(node_id, start, end))

    def get_subset(self, node_id, start, end):
        assert start >= 0
        assert end <= self.indexes(node_id)[-1]
        indexes = self.indexes(node_id)
        if len(indexes) == 2:
            return SimpleValuedIndexes(self.indexes, self.values)

        values = self.values(node_id)

        length = end-start
        indexes = indexes[1:-1]-start
        above_mask = indexes > 0
        above = np.nonzero(above_mask)[0]

        if not above.size:
            start_value = values[-1] if values.size else values[0]
        elif above[0] == 0:
            start_value = values[0]
        else:
            start_value = self.values[1:][above[0]-1]

        if not above.size:
            return SimpleValuedIndexes(np.array([0, length], dtype="int"),
                                       np.array([start_value]))

        inside_mask = np.logical_and(above_mask,
                                     indexes < length)
        subset_indexes = indexes[inside_mask]

        """
        if subset_indexes[0] != 0:
            subset_indexes = np.insert(subset_indexes, 0, 0)
        if subset_indexes[-1] != length:
            subset_indexes = np.append(subset_indexes, length)
        """

        subset_values = values[1:][inside_mask]
        subset_values = np.insert(subset_values, 0, start_value)

        return SimpleValuedIndexes(subset_indexes, subset_values)

    def get_max_value_between(self, node_id, start, end):
        assert node_id is not None
        assert start >= 0
        assert end <= self.indexes(node_id)[-1]

        max_val = 0
        for s, e, value in self.index_value_pairs(node_id):
            if s <= start and e > start:
                max_val = max(max_val, value)
            elif end > s and end <= e:
                max_val = max(max_val, value)

        return max_val

    def find_valued_areas(self, node, value):
        all_indexes = self.indexes(node)
        values = self.values(node)
        idxs = np.where(values == value)[0]
        starts = all_indexes[idxs]
        ends = all_indexes[idxs+1]
        areas = list(chain(*zip(starts, ends)))
        return areas

    def threshold_copy(self, cutoff):
        new = self.copy()
        new._values = new._values >= cutoff
        logging.info("Thresholding done.")
        return new

    def threshold(self, cutoff):
        if self._values is None:
            return  # Special case, empty pileup

        self._values = self._values >= cutoff

    def index_value_pairs(self, node):
        indexes = self.indexes(node)
        assert len(indexes) >= 2
        values = self.values(node)
        assert len(values) >= 1
        lines = list(zip(
            indexes[:-1],
            indexes[1:],
            values
            ))
        assert len(lines) >= 1
        return lines

    def get_bed_graph_lines(self):
        for node in self.nodes:
            lines = self.index_value_pairs(node)
            for line in lines:
                yield "%s\t%s\t%s\t%s\n" % (node, line[0], line[1], line[2])


class SparsePileup(Pileup):
    def __init__(self, graph):
        logging.info("Initing sparsepileup")
        self.graph = graph
        self.data = None

    def __eq__(self, other):
        return self.data == other.data

    def equals_old_sparse_pileup(self, old_pileup):
        # For tests to pass temporarily
        for node in self.data.nodes:
            indexes = self.data.indexes(node)
            other_indexes = old_pileup.data[node].all_idxs()
            if not np.all(indexes == other_indexes):
                return False

            values = self.data.values(node)
            other_values = old_pileup.data[node].all_values()
            if not np.allclose(values, other_values):
                return False

        return True

    @classmethod
    def from_base_value(cls, graph, base_value):
        pileup = cls(graph)
        pileup.data = SparsePileupData([], [], graph=graph, default_value=base_value)
        return pileup

    def sum(self):
        return self.data.sum()

    def mean(self):
        graph_size = sum([self.graph.node_size(b) for b in self.graph.blocks])
        mean = self.sum() / graph_size
        return mean

    def scale(self, scale):
        self.data.scale()

    def fill_small_wholes(self, max_size, write_holes_to_file=None, touched_nodes=None):
        cleaner = HolesCleaner(self, max_size, touched_nodes=touched_nodes)
        areas = cleaner.run()
        n_filled = 0

        hole_intervals = []

        for node_id in areas.areas:
            if touched_nodes is not None:
                if node_id not in touched_nodes:
                    continue

            starts = areas.get_starts(node_id)
            ends = areas.get_ends(node_id)
            for start, end in zip(starts, ends):
                logging.debug("Filling hole %s, %d, %d. Node size: %d" % (
                    node_id, start, end, self.graph.node_size(node_id)))

                if start == end:
                    logging.warning("Trying to fill hole of 0 length")
                    continue

                self.data.fill_existing_hole(node_id, start, True)

                n_filled += 1
                assert end - start <= max_size
                hole_intervals.append(Interval(start, end, [node_id]))

        logging.info(
            "Filled %d small holes (splitted into holes per node)" % n_filled)

        if write_holes_to_file is not None:
            intervals = IntervalCollection(hole_intervals)
            intervals.to_file(write_holes_to_file, text_file=True)

        self.sanitize()

    def sanitize(self):
        logging.info("Sanitizing sparse pileup")
        for node in self.data.nodes:
            self.data.sanitize_node(node)
        logging.info("Sanitize done")

    def find_valued_areas(self, value, only_check_nodes=None):
        if value:
            return {node_id: self.data.find_valued_areas(node_id, value)
                    for node_id in self.data.nodes}
        else:

            if only_check_nodes is not None:
                nodes = only_check_nodes
            else:
                nodes = self.graph.blocks

            return SparseAreasDict({node_id: self.data.find_valued_areas(node_id, value)
                                   for node_id in self.data.nodes
                                    }, graph=self.graph)

    def set_valued_intervals(self, node_id, valued_indexes):
        pass

    @classmethod
    def from_areas_collection(cls, graph, areas_list):
        starts_dict = defaultdict(list)
        ends_dict = defaultdict(list)
        for areas in areas_list:
            for rp in areas.areas:
                starts_dict[rp].extend(areas.get_starts(rp))
                ends_dict[rp].extend(areas.get_ends(rp))
        starts_dict = {rp: np.array(v) for rp, v in starts_dict.items()}
        ends_dict = {rp: np.array(v) for rp, v in ends_dict.items()}
        return cls.from_starts_and_ends(graph, starts_dict,
                                        ends_dict, dtype=int)


    @classmethod
    def init_sparsepileupdata_from_valued_ares(cls, graph, nodes, valued_areas):
        # Creates empty sparsepileupdata with correct size for each node
        data_lengths = []
        node_size = graph.blocks.node_size
        i = 0
        for rp in nodes:
            if i % 100000 == 0:
                logging.info("Initing sparse from valued areas for node %d" % i)
            i += 1

            length = node_size(rp)
            starts = valued_areas.get_starts_array(rp, node_size=length)
            ends = valued_areas.get_ends_array(rp, node_size=length)
            indexes, values = starts_and_ends_to_sparse_pileup(
                starts,
                ends)
            length = len(indexes) + 2
            data_lengths.append(length)

        pileup_data = SparsePileupData(nodes, data_lengths)
        del data_lengths
        return pileup_data

    @classmethod
    def from_valued_areas(cls, graph, valued_areas, touched_nodes = None):
        pileup = cls(graph)

        if touched_nodes is None:
            nodes = graph.blocks
        else:
            nodes = touched_nodes

        pileup_data = cls.init_sparsepileupdata_from_valued_ares(graph, nodes, valued_areas)

        i = 0
        # Fill pileup_data
        logging.info("N nodes to process: %d" % len(nodes))
        for rp in nodes:
            if i % 100000 == 0:
                logging.info("Creating sparse from valued areas for node %d" % i)
            i += 1

            length = graph.blocks[rp].length()
            starts = valued_areas.get_starts_array(rp, node_size=length)
            if len(starts) == 0:
                continue

            ends = valued_areas.get_ends_array(rp, node_size=length)
            if len(starts) == 0 and len(ends) == 0:
                continue

            indexes, values = starts_and_ends_to_sparse_pileup(
                starts,
                ends)

            if not indexes.size:
                #assert not valued_areas.has_anything_on_node(rp)
                continue

            start_value = 0
            #print("Orig indexes: %s" % indexes)
            if len(indexes) > 0:
                if indexes[-1] == length:
                    values = values[:-1]
                    indexes = indexes[:-1]

                if indexes[0] == 0:
                    start_value = values[0]
                    indexes = indexes[1:]
                    values = values[1:]

            #print("Indexes: %s" % indexes)
            #print("Length: %d" % length)
            pileup_data.set_indexes(rp, indexes)
            pileup_data.set_end_index(rp, length)
            pileup_data.set_start_value(rp, start_value)
            pileup_data.set_values(rp, values)
            #if len(indexes) > 0:
            #    print("Indexes/values:")
            #    print(pileup_data.indexes(rp))
            #    print(pileup_data.values(rp))

        pileup.data = pileup_data
        return pileup

    @classmethod
    def from_intervals(cls, graph, intervals):
        starts, ends = intervals_to_start_and_ends(graph, intervals)
        return cls.from_starts_and_ends(graph, starts, ends)

    @classmethod
    def from_starts_and_ends(cls, graph, starts, ends, dtype=bool):
        pileup = cls(graph)

        # Memory efficient and slower: First create empty structure with correct size, then fill
        node_ids = []
        lengths = []
        for rp in starts:
            assert rp in graph.blocks
            indexes, values = starts_and_ends_to_sparse_pileup(
                starts[rp], ends[rp])
            node_ids.append(rp)

            length = len(indexes) + 2
            if indexes[-1] == graph.blocks[rp].length():
                length -= 1
            if indexes[0] == 0:
                length -= 1
            lengths.append(length)

        data = SparsePileupData(node_ids, lengths, graph=graph)

        for rp in starts:
            assert rp in graph.blocks
            indexes, values = starts_and_ends_to_sparse_pileup(
                starts[rp], ends[rp])
            length = graph.blocks[rp].length()

            start_value = 0
            if len(indexes) > 0:
                if indexes[-1] == length:
                    values = values[:-1]
                    indexes = indexes[:-1]

                if indexes[0] == 0:
                    start_value = values[0]
                    indexes = indexes[1:]
                    values = values[1:]

            data.set_indexes(rp, indexes)
            data.set_end_index(rp, length)
            data.set_start_value(rp, start_value)
            data.set_values(rp, values)

        pileup.data = data

        logging.info("Number of elements in pileup: %d" % len(pileup.data.nodes))

        return pileup

    def to_subgraphs(self):
        # Returns a list of areas which each is a subgraph
        collection = SubgraphCollection.from_pileup(self.graph, self)
        return collection

    def __str__(self):
        out = "SparsepileupV2 \n"
        out += '\n'.join("  %d: %s, %s" % (rp, self.data.indexes(rp), self.data.values(rp)) \
                         for rp in self.data.nodes)
        return out
    __repr__ = __str__

    def threshold_copy(self, cutoff):
        new_pileup = self.__class__(self.graph)
        new_pileup.data = self.data.threshold_copy(cutoff)
        new_pileup.sanitize()
        return new_pileup

    def threshold(self, cutoff):
        self.data.threshold(cutoff)
        self.sanitize()

    def to_bed_file(self, filename):
        f = open(filename, "w")
        areas = self.find_valued_areas(True)
        for node_id, idxs in areas.items():
            for i in range(len(idxs)//2):
                interval = (node_id, idxs[2*i], idxs[2*i+1])
                f.write("%s\t%s\t%s\t.\t.\t.\n" % interval)
        f.close()
        return filename

    def remove_small_peaks(self, min_size):

        logging.info("Initing cleaner")
        cleaner = PeaksCleaner(self, min_size)
        logging.info("Running cleaner")
        areas = cleaner.run()
        logging.info("Removing emty areas")

        logging.info("Creating pileup using results from cleaner")
        pileup = self.from_areas_collection(self.graph, [areas])
        logging.info("Tresholding")
        pileup.threshold(0.5)
        return pileup

    def to_bed_graph(self, filename):
        f = open(filename, "w")
        i = 0
        for line in self.data.get_bed_graph_lines():
            f.write(line)
            i += 1
        f.close()
        self.filename = filename
        self.is_written = True
        return filename

    def set_sorted_interval_values(self, intervals, values):
        raise NotImplementedError()

    @classmethod
    def from_pickle(cls, file_name, graph):
        with open("%s" % file_name, "rb") as f:
            data = pickle.loads(f.read())
            #assert isinstance(SparsePileupData, cls)
            obj = cls(graph)
            obj.data = data
            return obj

    def to_pickle(self, file_name):
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self.data, f)

    @classmethod
    def create_from_old_sparsepileup(cls, old_pileup):
        # TMP method for tests to pass
        graph = old_pileup.graph
        nodes = []
        lengths = []

        for node in old_pileup.data:
            nodes.append(node)
            lengths.append(len(old_pileup.data[node].all_idxs()))

        pileupdata = SparsePileupData(nodes, lengths, graph=graph)
        for node in old_pileup.data:
            values = old_pileup.data[node].all_values()
            indexes = old_pileup.data[node].all_idxs()

            pileupdata.set_indexes(node, indexes[1:-1])
            pileupdata.set_end_index(node, indexes[-1])
            pileupdata.set_values(node, values[1:])
            pileupdata.set_start_value(node, values[0])

        pileup = SparsePileup(graph)
        pileup.data = pileupdata
        pileup.data.default_value = old_pileup

        return pileup