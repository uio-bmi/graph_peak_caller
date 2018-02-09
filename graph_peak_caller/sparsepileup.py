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


class ValuedIndexes(object):
    def __init__(self, indexes, values, start_value, length):

        if isinstance(indexes, list):
            indexes = np.array(indexes, dtype="int")

        if isinstance(values, list):
            values = np.array(values)


        assert type(indexes) == np.ndarray
        assert type(values) == np.ndarray
        # assert indexes.size == values.size
        self.length = length
        self.__values = None
        self.__indexes = None
        self.__start_value = 0
        self._all_values = np.zeros(0)
        self._all_indexes = np.zeros(0)

        self.values = values
        self.indexes = indexes
        self.start_value = start_value

        self._all_values[0] = start_value

        self.__tmp_starts = []
        self.__tmp_values = []
        self.__tmp_end = 0

    @property
    def indexes(self):
        #assert np.all(self.__indexes == self._all_indexes[1:-1])
        return self.__indexes

    @indexes.setter
    def indexes(self, new_indexes):
        self.__indexes = new_indexes
        if len(new_indexes) != len(self._all_indexes) - 2:
            self._all_indexes = np.zeros(len(new_indexes)+2, dtype=np.int)
        self._all_indexes[1:-1] = new_indexes
        self._all_indexes[-1] = self.length

    @property
    def values(self):
        return self.__values

    @values.setter
    def values(self, new_values):
        self.__values = new_values
        if len(new_values) != len(self._all_values) - 1 or new_values.ndim != self._all_values.ndim:
            if new_values.ndim == 2:
                self._all_values = np.zeros((len(new_values)+1, 2))
            else:
                self._all_values = np.zeros(len(new_values)+1)

            self._all_values[0] = self.start_value

        self._all_values[1:] = new_values

    @property
    def start_value(self):
        return self.__start_value

    @start_value.setter
    def start_value(self, value):
        self._all_values[0] = value
        self.__start_value = value

    def set_single_value(self, index, value):
        self.__values[index] = value
        self._all_values[index+1] = value

    def all_values(self):
        return self._all_values
        #return np.insert(self.values, 0, self.start_value)

    def all_idxs(self):
        return self._all_indexes
        #return np.append(np.insert(self.indexes, 0, 0), self.length)

    def __eq__(self, other):
        if not np.allclose(self.values, other.values):
            return False

        if not np.allclose(self.indexes, other.indexes):
            return False

        if isinstance(self.start_value, np.ndarray):
            if not np.allclose(self.start_value, other.start_value):
                return False
        else:
            if not np.isclose(self.start_value, other.start_value):
                return False
        return self.length == other.length

    def __str__(self):
        return "VI(indices: %s, values: %s, start: %s, length: %s)" % (
            self.indexes, self.values,
            self.start_value, self.length)

    __repr__ = __str__

    def sum(self):
        lengths = np.diff(self.all_idxs())
        return np.sum(lengths*self.all_values())

    def max_value(self):
        return np.max(self.all_values())

    def mean(self):
        return self.sum()/self.length

    def get_subset(self, start, end):
        assert start >= 0
        assert end <= self.length
        if not self.indexes.size:
            return ValuedIndexes(np.array([], dtype="int"),
                                 np.array([]), self.start_value,
                                 end-start)
        length = end-start
        indexes = self.indexes-start
        above_mask = indexes > 0
        above = np.nonzero(above_mask)[0]

        if not above.size:
            start_value = self.values[-1] if self.values.size else self.start_value
        elif above[0] == 0:
            start_value = self.start_value
        else:
            start_value = self.values[above[0]-1]
        if not above.size:
            return self.__class__(np.array([], dtype="int"),
                                  np.array([]),
                                  start_value,
                                  length)
        inside_mask = np.logical_and(above_mask,
                                     indexes < length)
        subset_indexes = indexes[inside_mask]
        subset_values = self.values[inside_mask]
        return self.__class__(subset_indexes,
                              subset_values, start_value,
                              length)

    def scale(self, scale):
        self.values = self.values*scale
        self.start_value = self.start_value*scale

    def tmp_set_interval_value(self, start, end, value):
        self.__tmp_starts.append(start)
        self.__tmp_values.append(value)
        self.__tmp_end = max(self.__tmp_end, end)

    def fix_tmp_values(self):
        assert self.__tmp_starts[0] == 0
        self.start_value = self.__tmp_values[0]
        self.length = self.__tmp_end
        if len(self.__tmp_starts) < 2:
            self.values = np.array([])
            self.indexes = np.array([], dtype="int")
        else:
            self.indexes = np.array(self.__tmp_starts[1:])
            self.values = np.array(self.__tmp_values[1:])
        self.__tmp_values = []
        self.__tmp_starts = []
        self.__tmp_end = 0

    def set_interval_value(self, start, end, value):
        # ONly used for filling holes. TODO rename
        if start == 0:
            self.start_value = value
        else:
            idx = np.nonzero(self.indexes == start)
            assert len(idx) == 1
            idx = idx[0]
            #self.values[idx] = value
            self.set_single_value(idx, value)

            assert end == self.length or np.any(self.indexes[idx[0]+1] == end)

    def set_interval_value_on_right_empty_area(self, start, end, value):
        # Requires that everything after the end is is empty
        if start == 0:
            self.start_value = value
        else:

            idx = np.where(start == self.indexes)[0]
            assert len(idx) <= 1
            if len(idx) == 1:
                self.values[idx] = value
            else:
                self.indexes = np.append(self.indexes, [start])
                self.values = np.append(self.values, [value])

        if end < self.length:
            self.indexes = np.append(self.indexes, [end])
            self.values = np.append(self.values, [0])

    def threshold_copy(self, cutoff):
        new_vi = self.__class__(np.copy(self.indexes),
                                self.values >= cutoff,
                                self.start_value >= cutoff,
                                self.length)
        new_vi.sanitize()
        return new_vi

    def threshold(self, cutoff):
        self.values = self.values >= cutoff
        self.start_value = self.start_value >= cutoff

    def sanitize(self):
        #diffs = np.diff(np.insert(self.values, 0, self.start_value))
        #assert np.all(self.all_values() == np.insert(self.values, 0, self.start_value))
        all_vals = self.all_values()

        diffs = np.diff(all_vals)
        changes = np.where(diffs != 0)[0]
        if len(self.values) == 0:
            return

        self.values = self.values[changes]
        self.indexes = self.indexes[changes]

    def sanitize_indices(self):
        indexes = self.all_idxs()
        values = self.all_values()
        diffs = np.where(np.diff(indexes) > 0)[0]

        self.indexes = indexes[diffs[1:]]
        self.values = values[diffs[1:]]
        self.start_value = values[diffs[0]]

    def find_valued_areas(self, value):
        all_indexes = self.all_idxs()
        values = self.all_values()
        idxs = np.where(values == value)[0]
        starts = all_indexes[idxs]
        ends = all_indexes[idxs+1]
        areas = list(chain(*zip(starts, ends)))
        return areas

    @classmethod
    def maximum(cls, vi_a, vi_b):
        a = vi_a.all_idxs()[:-1]*2
        b = vi_b.all_idxs()[:-1]*2+1
        all_idxs = np.concatenate([a, b])
        all_idxs.sort()
        vi_list = [vi_a, vi_b]
        values_list = []
        for i, vi in enumerate(vi_list):
            idxs = np.nonzero((all_idxs % 2) == i)[0]
            all_values = vi.all_values()
            value_diffs = np.diff(all_values)
            values = np.zeros(all_idxs.shape)
            values[idxs[1:]] = value_diffs
            values[idxs[0]] = all_values[0]
            values_list.append(values.cumsum())
        values = np.maximum(values_list[0], values_list[1])
        idxs = all_idxs // 2
        empty_ends = np.nonzero(np.diff(idxs) == 0)[0]
        max_values = np.maximum(values[empty_ends], values[empty_ends+1])
        values[empty_ends+1] = max_values
        values[empty_ends] = max_values
        obj = cls(idxs[1:], values[1:], values[0], vi_a.length)
        obj.sanitize()
        return obj

    @classmethod
    def combine(cls, vi_a, vi_b):
        a = vi_a.all_idxs()[:-1]*2
        b = vi_b.all_idxs()[:-1]*2+1
        all_idxs = np.concatenate([a, b])
        all_idxs.sort()
        vi_list = [vi_a, vi_b]
        values_list = []
        for i, vi in enumerate(vi_list):
            idxs = np.nonzero((all_idxs % 2) == i)[0]
            all_values = vi.all_values()
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
        obj = cls(idxs[1:], np.transpose(values[:, 1:]),
                  values[:, 0], vi_a.length)
        return obj

    def trunctate(self, min_value):
        self.values = np.maximum(self.values, min_value)
        self.start_value = max(min_value, self.start_value)
        self.sanitize()

    @classmethod
    def empty(cls, length):
        return cls(np.array([], dtype="int"),
                   np.array([]),
                   0,
                   length)

    def __iter__(self):
        return zip(
            chain([0], self.indexes),
            chain(self.indexes, [self.length]),
            chain([self.start_value], self.values)
            )


class BinaryIndexes(object):
    def __init__(self, starts, ends, length):
        self.starts = starts
        self.ends = ends
        self.length = length

    def add_interval(self, start, end):
        self.starts.append(start)
        self.ends.append(end)
        self.starts.sort


class SparsePileupData(dict):

    def __init__(self, *args, **kwargs):
        self.graph = kwargs["graph"]
        self.min_value = 0
        del kwargs["graph"]
        super(SparsePileupData, self).__init__(*args, **kwargs)
        self.default_value = 0

    def __getitem__(self, item):
        if item not in self:
            value = ValuedIndexes.empty(self.graph.node_size(item))
            value.start_value = max(self.min_value, self.default_value)
            self.__setitem__(item, value)

        return super(SparsePileupData, self).__getitem__(item)

    def __setitem__(self, key, value):
        super(SparsePileupData, self).__setitem__(key, value)

    def trunctate(self, min_value):
        for key in self.keys():
            self.__getitem__(key).trunctate()
        self.min_value = max(self.min_value, min_value)

    def values(self):
        return (self.__getitem__(key) for key in self.graph.blocks)

    # def items(self):
    #    return ((key, self.__getitem__(key)) for key in self.graph.blocks)

    def all_items(self):
        # Also gives rps with zero
        return ((key, self.__getitem__(key)) for key in self.graph.blocks)


class NpSparseAreasDict(dict):
    def __init__(self, diffs, changes):
        self._diffs = diffs
        self._changes = changes

    def __getitem__(self, item):
        pass

class SparseAreasDict(dict):
    """
    Used as return object in SparsePileup.get_valued_areas() to speedup.
    Assumes value in empty nodes, i.e. nodes not in dict are completely filled
    """
    def __init__(self, *args, **kwargs):
        self.graph = kwargs["graph"]
        self.min_value = 0
        del kwargs["graph"]
        super(SparseAreasDict, self).__init__(*args, **kwargs)

    def __getitem__(self, item):
        if item not in self:
            value = [0, self.graph.node_size(item)]
            self.__setitem__(item, value)

        return super(SparseAreasDict, self).__getitem__(item)

    def values(self):
        raise NotImplementedError()

    def items(self):
        return ((key, self.__getitem__(key)) for key in self.graph.blocks)


class SparsePileup(Pileup):
    def __init__(self, graph):
        logging.info("Initing sparsepileup")
        self.graph = graph
        self.data = SparsePileupData(graph=self.graph)

        #self.data = {rp: ValuedIndexes.empty(graph.node_size(rp))
        #             for rp in self.graph.blocks}
        logging.info("Sparsepileup inited")
        #self.graph.assert_correct_edge_dicts()

    def __eq__(self, other):
        for node_id, vi in other.data.all_items():
            if self.data[node_id] != vi:
                return False
        return True

    @classmethod
    def from_base_value(cls, graph, base_value):
        pileup = cls(graph)
        pileup.update_max_value(base_value)
        return pileup

    def sum(self):
        return np.sum([values.sum() for node, values in self.data.items()])

    def mean(self):
        graph_size = sum([self.graph.node_size(b) for b in self.graph.blocks])
        mean = self.sum() / graph_size
        return mean

    def scale(self, scale):
        [vi.scale(scale) for vi in self.data.values()]

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
                self.data[node_id].set_interval_value(start, end, True)
                logging.debug("Filling hole %s, %d, %d" % (
                    node_id, start, end))
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
        for node_id, valued_indexes in self.data.items():
            valued_indexes.sanitize()
            assert node_id in self.graph.blocks

    def find_valued_areas(self, value, only_check_nodes=None):
        if value:
            return {node_id: self.data[node_id].find_valued_areas(value)
                    for node_id, valued_indexes in self.data.items()}
        else:

            if only_check_nodes is not None:
                nodes = only_check_nodes
            else:
                nodes = self.graph.blocks

            return SparseAreasDict({node_id: self.data[node_id].find_valued_areas(value)
                                   for node_id in self.data
                                    }, graph=self.graph)

            """
            return {node_id: (self.data[node_id].find_valued_areas(value)
                              if node_id in self.data
                              else [0, self.graph.node_size(node_id)])
                    for node_id in nodes}
            """

    def set_valued_intervals(self, node_id, valued_indexes):
        assert node_id in self.graph.blocks
        self.data[node_id] = valued_indexes

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
            if indexes[0] == 0:
                start_value = values[0]
                indexes = indexes[1:]
                values = values[1:]

            if len(indexes) > 0:
                if indexes[-1] == length:
                    indexes = indexes[:-1]
                    values = values[:-1]

            pileup.data[rp] = ValuedIndexes(
                indexes, values, start_value, length)

        return pileup

    @classmethod
    def from_starts_and_ends(cls, graph, starts, ends, dtype=bool):
        pileup = cls(graph)
        for rp in starts:
            assert rp in graph.blocks
            indexes, values = starts_and_ends_to_sparse_pileup(
                starts[rp], ends[rp])
            start_value = dtype(False)
            length = graph.blocks[rp].length()
            if indexes[0] == 0:
                start_value = values[0]
                indexes = indexes[1:]
                values = values[1:]

            if len(indexes) > 0:
                if indexes[-1] == length:
                    indexes = indexes[:-1]
                    values = values[:-1]

            pileup.data[rp] = ValuedIndexes(
                indexes, values, start_value, length)

        logging.info("Number of elements in pileup: %d" % len(pileup.data))
        return pileup

    @classmethod
    def from_intervals(cls, graph, intervals):
        starts, ends = intervals_to_start_and_ends(graph, intervals)
        return cls.from_starts_and_ends(graph, starts, ends)

    def to_subgraphs(self):
        # Returns a list of areas which each is a subgraph
        collection = SubgraphCollection.from_pileup(self.graph, self)
        return collection

    def __str__(self):
        return "\n".join(
            "%s: %s, %s, %s" % (node_id, valued_indexes.indexes, valued_indexes.values, valued_indexes.start_value)
            for node_id, valued_indexes in self.data.items())

    __repr__ = __str__

    def threshold_copy(self, cutoff):
        new_data = {node_id: vi.threshold_copy(cutoff)
                    for node_id, vi in self.data.items()}

        pileup = self.__class__(self.graph)
        pileup.data = new_data
        return pileup

    def threshold(self, cutoff):
        n = 0
        logging.info("Thresholding. Number of nodes to threshold: %d" % len(self.data))
        for node, valued_indexes in self.data.items():
            if n % 25000 == 0:
                logging.info("Thresholding node %d" % n)
            n += 1
            valued_indexes.threshold(cutoff)
            valued_indexes.sanitize()

    def add_intervals(self, intervals):
        [self.add_interval(interval)
         for interval in intervals]

    def add_interval(self, interval):
        assert all(region_path in self.graph.blocks for
                   region_path in interval.region_paths)
        for i, region_path in enumerate(interval.region_paths):
            start = 0
            end = self.graph.node_size(region_path)
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths)-1:
                end = interval.end_position.offset
            self.data[region_path].add_interval(
                start, end)

    def add_area(self, region_path, start, end, value):
        """
        NAIVE method to add single area to valued_inexes
        Assumes it comes from a complete set of areas
        """
        vi = self.data[region_path]
        vi.tmp_set_interval_value(start, end, value)

    def fix_tmp_values(self):
        logging.info("Fixing tmp values")
        for vi in self.data.values():
            vi.fix_tmp_values()

    @classmethod
    def from_bed_graph(cls, graph, filename):
        pileup = super().from_bed_graph(graph, filename)
        assert isinstance(pileup, cls)
        pileup.fix_tmp_values()
        return pileup

    def add_areas(self, areas):
        for area, intervals in areas.items():
            self.data[area].add_interval(
                intervals[0], intervals[1])

    def set_areas_value(self, areas, value):
        for area, intervals in areas.items():
            for i in range(len(intervals)//2):
                self.data[area].set_interval_value(
                    intervals[i], intervals[i+1], value)

    def set_interval_value(self, interval, value):
        assert all(region_path in self.graph.blocks for
                   region_path in interval.region_paths)
        for i, region_path in enumerate(interval.region_paths):
            start = 0
            end = self.graph.node_size(region_path)
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths)-1:
                end = interval.end_position.offset
            self.data[region_path].set_interval_value(
                start, end, value)

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

    def set_to_false(self):
        for valued_indexes in self.data.values():
            valued_indexes.start_value = False
            valued_indexes.indexes = np.array([], dtype="int")
            valued_indexes.values = np.array([], dtype="bool")

    def remove_small_peaks(self, min_size):
        """
        areas = self.find_valued_areas(True)
        intervals = self.areas_to_intervals(areas, include_partial_stubs=False)
        large_intervals = [interval for interval in intervals
                           if interval.length() >= min_size]
        for i in [interval for interval in intervals
                  if interval.length() < min_size]:
            print(i)
        return self.from_intervals(self.graph, large_intervals)
        """
        logging.info("Initing cleaner")
        cleaner = PeaksCleaner(self, min_size)
        logging.info("Running cleaner")
        areas = cleaner.run()
        logging.info("Removing emty areas")
        #for node, startends in areas.areas.items():
        #    if len(startends) == 0:
        #        del areas[node]

        logging.info("Creating pileup using results from cleaner")
        pileup = self.from_areas_collection(self.graph, [areas])
        logging.info("Tresholding")
        pileup.threshold(0.5)
        return pileup

    def update_max(self, other):
        for key, valued_indexes in self.data.items():
            self.data[key] = ValuedIndexes.maximum(
                valued_indexes, other.data[key])

    def update_max_value(self, min_value):

        if isinstance(self.data, SparsePileupData):
            self.data.trunctate(min_value)
        else:
            for valued_indexes in self.data.values():
                valued_indexes.trunctate(min_value)

    def to_bed_graph(self, filename):
        f = open(filename, "w")
        for node_id, valued_indexes in self.data.items():
            for t in valued_indexes:
                f.write("%s\t%s\t%s\t%s\n" % ((node_id,) + t))
        f.close()
        self.filename = filename
        self.is_written = True
        return filename

    def set_sorted_interval_values(self, intervals, values):
        # Requires intervals to be sorted within nodes
        for j, interval in enumerate(intervals):
            logging.info("Setting sorted interval %d/%d" % (j, len(intervals)))
            for i, rp in enumerate(interval.region_paths):
                start = 0
                if i == 0:
                    start = interval.start_position.offset
                end = self.graph.node_size(rp)
                if i + 1 == len(interval.region_paths):
                    end = interval.end_position.offset
                self.data[rp].set_interval_value_on_right_empty_area(start, end, values[j])

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
        logging.info("N nodes in control: %d" % len(control.data))
        sparse_pileup = cls(sample.graph)
        sparse_pileup.data = {
            node_id: ValuedIndexes.combine(
                control.data[node_id], sample.data[node_id])
            for node_id in control.data}
        return sparse_pileup

    @classmethod
    def from_control_and_sample(cls, control, sample):
        sparse_pileup = cls(control.graph)
        for node_id in sparse_pileup.graph.blocks:
            control_counts = control.get_count_array(node_id)
            sample_counts = sample.get_count_array(node_id)
            control_diffs = control_counts[1:]-control_counts[:-1]
            sample_diffs = sample_counts[1:]-sample_counts[:-1]
            changes = np.where(np.logical_or(control_diffs != 0,
                                             sample_diffs))[0] + 1
            vals = np.column_stack([control_counts[changes],
                                    np.floor(sample_counts[changes])])
            init_value = np.array([control_counts[0], sample_counts[0]])
            valued_indexes = ValuedIndexes(
                 changes, vals, init_value, control_counts.shape[0])
            assert not np.any(vals[:, 0] == 0)
            sparse_pileup.set_valued_intervals(node_id, valued_indexes)
        return sparse_pileup


def intervals_to_start_and_ends(graph, intervals):
    # Returns two dict on rp => positions (start/ends)
    starts = defaultdict(list)
    ends = defaultdict(list)

    for interval in intervals:
        for i, region_path in enumerate(interval.region_paths):
            start = 0
            end = graph.node_size(region_path)
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths)-1:
                end = interval.end_position.offset

            if region_path < 0:
                new_start = graph.node_size(region_path) - end
                new_end = graph.node_size(region_path) - start
                start = new_start
                end = new_end

                region_path = -region_path

            starts[region_path].append(start)
            ends[region_path].append(end)

    for rp in starts:
        starts[rp] = np.array(starts[rp])

    for rp in ends:
        ends[rp] = np.array(ends[rp])

    return starts, ends


def filter_pileup_duplicated_position(positions, values):
    equal_to_previous = (positions[:-1] == positions[1:])
    true_array = np.ones(len(positions))
    true_array[np.where(equal_to_previous)] = False
    return positions[np.where(true_array)], values[np.where(true_array)]


def starts_and_ends_to_sparse_pileup(starts, ends):
    indices, values = filter_pileup_duplicated_position(
        *DiscreteEventSorter([ends, starts]).pileup())
    return indices, values

