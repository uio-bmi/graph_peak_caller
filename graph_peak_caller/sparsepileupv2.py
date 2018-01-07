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
from .sparsepileup import SparseAreasDict, starts_and_ends_to_sparse_pileup


class SparsePileupData:
    def __init__(self, node_ids, lengths, ndim=1):

        self._node_indexes = None  # Mapping from node id to index in indexes/values
        self._indexes = None
        self._values = None
        self._lengths = None
        self.min_node = None
        self.min_value = 0

        self._create_empty(node_ids, lengths, ndim)

    def _create_empty(self, node_ids, lengths, ndim):

        self.nodes = set(node_ids)

        sorted_nodes = sorted(node_ids)
        self.min_node = sorted_nodes[0]
        max_node = sorted_nodes[-1]
        span = max_node - self.min_node + 1

        n_elements = sum(lengths)
        self._indexes = np.zeros(n_elements, dtype=np.uint32)
        if ndim > 1:
            self._values = np.zeros((n_elements, ndim))
        else:
            self._values = np.zeros(n_elements)
        self._node_indexes = np.zeros(span, dtype=np.uint32)
        self._lengths = np.zeros(span, dtype=np.uint32)

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
        new = SparsePileupData([], [])
        new._indexes = self.indexes.copy()
        new._values = self.values.copy()
        new._node_indexes= self._node_indexes.copy()
        new._lengths = self._lengths.copy()
        new.min_node = self.min_node.copy()

    def set_min_value(self, value):
        self.min_value = value

    def set_values(self, node_id, values):
        # First index is always start value
        # This method only sets everything else than start
        assert node_id in self.nodes, "Only allowed to set for already initiated nodes"
        index = node_id - self.min_node
        start = self._node_indexes[index] + 1
        end = start + len(values)
        self._values[start:end] = values

    def set_indexes(self, node_id, indexes):
        # First index is always 0. Last index is always length
        # This method sets everything else than first index
        assert node_id in self.nodes, "Only allowed to set for already initiated nodes"
        index = node_id - self.min_node
        start = self._node_indexes[index] + 1
        end = start + len(indexes)
        assert len(indexes) <= self._lengths[index]
        print(start)
        print(end)
        print(type(start))
        self._indexes[start:end] = indexes

        #print(self._values)
        #print(self._indexes)
        #print(self._lengths)
        #print(self._node_indexes)

        assert self._indexes[start-1] == 0, "Start index should be 0, is %d" % self._indexes[start-1]

        assert np.all(self.indexes(node_id)[1:-1] == indexes), "%s != %s" % (self.indexes(node_id)[1:-1], indexes)

    def set_end_index(self, node_id, i):
        # Used to set "length" of node. End should always be length
        index = node_id - self.min_node
        start = self._node_indexes[index]
        end = start + self._lengths[index]

        for idx in self.indexes(node_id):
            assert i > idx, "Trying to set end idx %d which is <= already existing index %d. All indexes: %s" % (i, idx, self.indexes(node_id))

        assert self._indexes[end-1] == 0, "End index already has something"
        self._indexes[end-1] = i
        assert self._indexes[end-1] == i, "%d != %d" % (self._indexes[end-1], i)

        assert self.indexes(node_id)[-1] == i, "End index %d is not %d. All indexes: %s" % \
            (i, self.indexes(node_id)[-1], self.indexes(node_id))


    def set_start_value(self, node_id, val):
        index = node_id - self.min_node
        start = self._node_indexes[index]
        self._values[start] = val

    def indexes(self, node_id):
        index = node_id - self.min_node

        if index < 0 or index >= len(self._node_indexes):
            return np.array([])

        length = self._lengths[index]
        assert length >= 2
        start = self._node_indexes[index]
        end = start + length
        return self._indexes[start:end]

    def values(self, node_id):
        index = node_id - self.min_node
        if index < 0 or index >= len(self._node_indexes):
            return np.array([])

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
        old_values = self.values(node)
        old_values[idx] = value
        self.set_values(node, old_values)

    def sanitize_node(self, node_id):
        indexes = self.indexes(node_id)
        values = self.values(node_id)
        diffs = np.where(np.diff(indexes) > 0)[0]

        self.set_indexes(indexes[diffs])
        self.set_values(indexes[diffs[1:]])
        self._lengths[node_id - self.min_node] = len(diffs)

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
        return new

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

    @classmethod
    def combine_valued_indexes(cls, indexes1, values1, indexes2, values2):

        a = indexes1[:-1]*2
        b = indexes2[:-1]*2+1
        all_idxs = np.concatenate([a, b])
        all_idxs.sort()

        values_list = []
        for i, vi in enumerate((values1, values2)):
            idxs = np.nonzero((all_idxs % 2) == i)[0]
            all_values = vi
            value_diffs = np.diff(all_values)
            values = np.zeros(all_idxs.shape)
            values[idxs[1:]] = value_diffs
            values[0] = all_values[0]
            values_list.append(values.cumsum())

        values = np.array([values_list[0], values_list[1]])
        idxs = all_idxs // 2
        unique_idxs = np.append(np.nonzero(np.diff(idxs))[0], len(idxs)-1)
        idxs = idxs[unique_idxs]
        values = values[:, unique_idxs]

        return (idxs, np.transpose(values))

        #obj = cls(idxs[1:], np.transpose(values[:, 1:]),
        #          values[:, 0], vi_a.length)
        #return obj

class SparsePileup(Pileup):
    def __init__(self, graph):
        logging.info("Initing sparsepileup")
        self.graph = graph
        self.data = None

    def __eq__(self, other):
        pass

    @classmethod
    def from_base_value(cls, graph, base_value):
        pileup = cls(graph)
        pileup.data.set_min_value(base_value)
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
                logging.debug("Filling hole %s, %d, %d" % (
                    node_id, start, end))

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
        for node in self.data.nodes:
            self.data.sanitize_node(node)

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
    def from_valued_areas(cls, graph, valued_areas, touched_nodes = None):
        pileup = cls(graph)
        i = 0

        if touched_nodes is None:
            nodes = graph.blocks
        else:
            nodes = touched_nodes

        data_nodes = []
        data_lengths = []

        # First create empty data
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

            length = len(indexes) + 2
            if indexes[-1] == graph.blocks[rp].length():
                length -= 1
            if indexes[0] == 0:
                length -= 1

            assert length >= 2

            data_nodes.append(rp)
            data_lengths.append(length)

        pileup_data = SparsePileupData(data_nodes, data_lengths)

        # Fill pileup_data
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

        data = SparsePileupData(node_ids, lengths)

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
        new_pileup.data = self.data.treshold_copy(cutoff)
        return new_pileup

    def add_area(self, region_path, start, end, value):
        """
        NAIVE method to add single area to valued_inexes
        Assumes it comes from a complete set of areas
        """
        raise NotImplementedError()

    def fix_tmp_values(self):
        raise NotImplementedError()

    @classmethod
    def from_bed_graph(cls, graph, filename):
        pileup = super().from_bed_graph(graph, filename)
        assert isinstance(pileup, cls)
        pileup.fix_tmp_values()
        return pileup

    @classmethod
    def from_bed_file(cls, graph, filename):
        f = open(filename, "r")
        starts = defaultdict(list)
        ends = defaultdict(list)
        for line in f:
            if line.startswith("track"):
                continue

            data = line.split()
            block_id = int(data[0])
            starts[block_id].append(int(data[1]))
            ends[block_id].append(int(data[2]))
        starts = {block_id: np.array(start_list) for block_id, start_list in starts.items() if start_list}
        ends = {block_id: np.array(ends_list) for block_id, ends_list in ends.items() if ends_list}
        return cls.from_starts_and_ends(graph, starts, ends)

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


class SparseControlSample(SparsePileup):
    def get_p_dict(self):
        p_value_dict = defaultdict(dict)
        count_dict = defaultdict(int)
        baseEtoTen = np.log(10)
        for node_id, valued_indexes in self.data.items():
            start = 0
            val = valued_indexes.start_value
            for start, end, val in valued_indexes:
                if val[1] not in p_value_dict[val[0]]:
                    log_p_val = poisson.logsf(val[1], val[0])
                    p_value_dict[val[0]][val[1]] = -log_p_val/baseEtoTen
                p = p_value_dict[val[0]][val[1]]
                count_dict[p] += end-start

        self.p_value_dict = p_value_dict
        self.count_dict = count_dict

    def get_p_to_q_values(self):
        p_value_counts = self.count_dict
        p_to_q_values = {}
        sorted_p_values = sorted(p_value_counts.keys(), reverse=True)
        rank = 1
        # logN = np.log10(self.graph.get_size())
        logN = np.log10(sum(p_value_counts.values()))
        pre_q = None
        for p_value in sorted_p_values:
            value_count = p_value_counts[p_value]
            q_value = p_value + (np.log10(rank) - logN)
            if np.isclose(p_value, 2.1326711212014025):
                print(q_value, logN, rank, np.log10(rank))
            if rank == 1:
                q_value = max(0.0, q_value)
            else:
                q_value = max(0.0, min(pre_q, q_value))
            p_to_q_values[p_value] = q_value
            pre_q = q_value
            rank += value_count
        self.p_to_q_values = p_to_q_values

    def get_q_values(self):
        def translation(x):
            return self.p_to_q_values[
                self.p_value_dict[x[0]][x[1]]]
        # f = np.vectorize(translation)
        for valued_indexes in self.data.values():

            valued_indexes.start_value = translation(
                valued_indexes.start_value)

            if valued_indexes.values.size:
                new_values = np.apply_along_axis(
                    translation, 1, valued_indexes.values)
                valued_indexes.values = new_values

            valued_indexes.sanitize()

    def get_scores(self):
        self.get_p_dict()
        self.get_p_to_q_values()
        self.get_q_values()

    @classmethod
    def from_sparse_control_and_sample(cls, control, sample):
        logging.info("Creating pileup by combining control and sample")

        # Create new empty sparse pileup with enough space for each node (sum of lengths)
        all_nodes = list(control.data.nodes.union(sample.data.nodes))
        lengths = []
        for node in all_nodes:
            all_indexes = np.append(control.data.indexes(node), sample.data.indexes(node))
            unique = np.unique(all_indexes)
            print("UNique: %s" % unique)
            lengths.append(len(unique))
            print("Node %d has length %d" % (node, len(unique)))

        pileupdata = SparsePileupData(all_nodes, lengths, ndim=2)

        sparse_pileup = cls(sample.graph)
        for node in all_nodes:
            indexes, value_pairs = SparsePileupData.combine_valued_indexes(
                control.data.indexes(node),
                control.data.values(node),
                sample.data.indexes(node),
                sample.data.values(node)
            )
            print("Combined values:")
            print(value_pairs)
            print("Combined indices:")
            print(indexes)

            pileupdata.set_indexes(node, indexes[1:])
            pileupdata.set_end_index(node, sample.graph.node_size(node))
            pileupdata.set_values(node, value_pairs[1:, :])
            pileupdata.set_start_value(node, value_pairs[0, :])

        sparse_pileup.data = pileupdata

        return sparse_pileup



