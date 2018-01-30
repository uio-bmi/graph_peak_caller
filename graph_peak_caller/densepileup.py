import logging
import numpy as np
from scipy.stats import poisson
from collections import defaultdict
from .pileup import Pileup
from .pileupcleaner2 import PeaksCleaner, HolesCleaner
from .subgraphcollection import SubgraphCollection
from offsetbasedgraph import Interval, IntervalCollection, BlockArray
import pickle
from .sparsepileup import SparseAreasDict, intervals_to_start_and_ends
# from memory_profiler import profile
from .sparsepileupv2 import RpScore


class DensePileupData:

    def __init__(self, graph, ndim=1, base_value=0, dtype=None):
        self._values = None
        self._node_indexes = None
        self._graph = graph
        self.min_node = None
        self._touched_nodes = set()
        self.ndim = ndim
        self.dtype = dtype

        self._create_empty(ndim, base_value=base_value)

    def _create_empty(self, ndim=1, base_value=0):
        logging.info("Sorting nodes")
        self._nodes = sorted(self._graph.blocks.keys())
        sorted_nodes = self._nodes
        self.min_node = sorted_nodes[0]
        max_node = sorted_nodes[-1]
        span = max_node - self.min_node + 1
        logging.info("Counting basepairs")
        n_elements = self._graph.number_of_basepairs()

        if self.dtype is not None:
            self._values = np.zeros(n_elements)
        else:
            self._values = np.zeros(n_elements)

        if base_value > 0:
            self._values += base_value

        if isinstance(self._graph.blocks, BlockArray):
            # Quicker way to make node_indexes array
            logging.info("Using quick way to init densepileup (using cumsum on np block array)")
            self._node_indexes = np.cumsum(self._graph.blocks._array, dtype=np.uint32)
            logging.info("Node indexes created...")
        else:
            self._node_indexes = np.zeros(span, dtype=np.uint32)
            offset = 0
            for i, node in enumerate(self._nodes):
                index = node - self.min_node
                self._node_indexes[index] = offset
                offset += self.node_size(node)

        logging.info("Dense pileup inited")

    def sum(self):
        return np.sum(self._values)

    def sanitize_node(self):
        return

    def node_size(self, node):
        return self._graph.node_size(node)

    def values(self, node):
        index = node - self.min_node
        start = self._node_indexes[index]
        end = start + self.node_size(node)
        return self._values[start:end]

    def values_in_range(self, node, start, end):
        index = node - self.min_node
        array_start = self._node_indexes[index] + start
        array_end = self._node_indexes[index] + end
        return self._values[array_start:array_end]

    def node_range_to_value_indexes(self, node, start, end):

        if node < 0:
            node = -node
            new_start = self.node_size(node) - end
            end = self.node_size(node) - start
            start = new_start

        node_size = self.node_size(node)
        start = min(node_size, max(start, 0))
        end = max(0, min(node_size, end))

        index = node - self.min_node
        array_start = self._node_indexes[index] + start
        array_end = self._node_indexes[index] + end
        return array_start, array_end

    def set_values(self, node, start, end, value):
        index = node - self.min_node
        array_start = self._node_indexes[index] + start
        array_end = self._node_indexes[index] + end
        self._values[array_start:array_end] = value
        self._touched_nodes.add(node)

    def add_value(self, node, start, end, value):
        index = node - self.min_node
        array_start = self._node_indexes[index] + start
        array_end = self._node_indexes[index] + end
        self._values[array_start:array_end] += value
        self._touched_nodes.add(node)

    def add_value_to_full_node(self, node, value):
        index = node - self.min_node
        array_start = self._node_indexes[index]
        array_end = array_start + self.node_size(node)
        self._values[array_start:array_end] += value
        self._touched_nodes.add(node)

    def set_full_node_value(self, node, value):
        index = node - self.min_node
        array_start = self._node_indexes[index]
        array_end = array_start + self.node_size(node)
        self._values[array_start:array_end] = value
        self._touched_nodes.add(node)

    def get_subset_max_value(self, node_id, start, end):
        return np.max(self.values(node_id)[start:end])

    def get_subset_sum(self, node_id, start, end):
        return np.sum(self.values(node_id)[start:end])

    def score(self, node_id, start, end):
        return RpScore(self.get_subset_max_value(node_id, start, end), \
               self.get_subset_sum(node_id, start, end))

    def scale(self, factor):
        self._values *= factor

    def set_new_values(self, values):
        self._values = values

    def fill_existing_hole(self, node, start, end, value):
        assert np.all(self.values_in_range(node, start, end) == 0)
        self.set_values(node, start, end, value)

    def get_sparse_indexes_and_values(self, node):
        values = self.values(node)
        diffs = np.ediff1d(values, to_begin=np.array([1]))
        indexes = np.where(diffs != 0)
        values = values[indexes]
        indexes = np.append(indexes, self.node_size(node))
        return indexes, values

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

    def find_valued_areas(self, node, value, changes=None):
        # Return list[start, end, start2, end2,...] having this value inside
        index = node - self.min_node
        start = self._node_indexes[index]
        end = start + self.node_size(node)
        changes = np.nonzero(changes[start:end])[0]+1
        # is_value = values == value
        # changes = np.nonzero(np.ediff1d(is_value))[0]+1
        if self._values[start] == value:
            if self._values[end-1] == value:
                return [0]+list(changes)+[end-start]
            return [0]+list(changes)
        if self._values[end-1]:
            return list(changes)+[end-start]
        return list(changes)

    def nodes(self):
        return self._touched_nodes

    def copy(self):
        new = DensePileupData(self._graph)
        new._values = self._values
        new._touched_nodes = self._touched_nodes
        return new

    def threshold_copy(self, cutoff):
        new = self.copy()
        new._values = new._values >= cutoff
        logging.info("Thresholding done.")
        return new

    def threshold(self, cutoff):
        self._values = self._values >= cutoff

    def index_value_pairs(self, node):
        indexes, values = self.get_sparse_indexes_and_values(node)
        assert len(indexes) >= 2
        assert len(values) >= 1
        lines = list(zip(
            indexes[:-1],
            indexes[1:],
            values
            ))
        assert len(lines) >= 1
        return lines

    def get_bed_graph_lines(self):
        for node in self.nodes():
            lines = self.index_value_pairs(node)
            for line in lines:
                yield "%s\t%s\t%s\t%s\n" % (node, line[0], line[1], line[2])

    def __len__(self):
        return len(self.nodes())

    def __eq__(self, other):
        for node in self.nodes():
            indexes, values = self.get_sparse_indexes_and_values(node)
            #print("   Checking node %d" % node)
            other_indexes, other_values = other.get_sparse_indexes_and_values(node)
            #print("   Values: %s / %s" % (values, other_values))
            if not np.all(indexes == other_indexes):
                print("Indices %s != %s" % (indexes, other_indexes))
                return False


            if np.any(np.abs(values -  other_values) > 1e-5):
                print("Values %s != %s" % (values, other_values))
                print()
                return False

        return True

    def get_interval_values(self, interval):

        values = np.zeros(interval.length())
        offset = 0
        """
        find_reversed = False
        use_interval = interval

        if np.all([rp < 0 for rp in interval.region_paths]):
            # Reverse before finding
            find_reversed = True
            use_interval = interval.get_reverse()
        else:
            assert np.all([rp > 0 for rp in interval.region_paths]), \
                "This method only supports intervals with single rp direction"
        """
        for i, rp in enumerate(interval.region_paths):
            assert rp > 0, "Currently only implemented for forward directed intervals"
            start = 0
            end = self._graph.node_size(rp)
            if i == 0:
                start = interval.start_position.offset
            if i == len(interval.region_paths) - 1:
                end = interval.end_position.offset

            values_in_rp = self.values_in_range(rp, start, end)
            values[offset:offset + (end - start)] = values_in_rp

            offset += end-start

        return values

    def value_indexes_to_nodes_and_offsets(self, indexes):
        """
        Takes indexes referring to positions in self._values
        and returns list of (node, offset).
        Indexes must be sorted
        """

        positions = []
        # Map indexes to node and offset
        length_offset = 0
        nodes = self._nodes
        current_node_index = 0
        i = 0
        while i < len(indexes):
            index = indexes[i]
            current_node = nodes[current_node_index]
            node_size = self.node_size(current_node)
            next_node_start = length_offset + node_size

            #print("Checking index %d. Current node: %d" % (index, current_node))

            if index >= next_node_start:
                current_node_index += 1
                length_offset += node_size
                continue

            if index < next_node_start:
                positions.append((current_node, index - length_offset))

            i += 1
        return positions


class DensePileup(Pileup):
    def __init__(self, graph, ndim=1, base_value=0, dtype=None):
        logging.info("Initing dense pileup")
        self.graph = graph
        self.data = DensePileupData(graph, ndim=ndim, base_value=base_value, dtype=dtype)

    def add_areas(self, areas):
        for area, intervals in areas.items():
            self.data.values(area)[intervals[0]:intervals[1]] += 1
        self.data._touched_nodes.add(area)

    @classmethod
    def from_base_value(cls, graph, base_value):
        pileup = cls(graph, base_value=base_value)
        return pileup

    def __str__(self):
        out = "Densepileup \n"
        for node in self.data._touched_nodes:
            #out += "  Node %d: %s, %s\n" % (node, self.data.values(node), self.data.get_sparse_indexes_and_values(node))
            out += "  Node %d: %s\n" % (node, self.data.get_sparse_indexes_and_values(node))

        return out

    def __repr__(self):
        return self.__str__

    def __eq__(self, other):
        return self.data == other.data

    def sum(self):
        return self.data.sum()

    def mean(self):
        graph_size = sum([self.graph.node_size(b) for b in self.graph.blocks])
        mean = self.sum() / graph_size
        return mean

    def scale(self, factor):
        self.data.scale(factor)

    def fill_small_wholes_on_dag(self, max_size):
        from .dagholecleaner import DagHoleCleaner
        logging.info("Cleaning holes using Dag Hole Cleaner")
        cleaner = DagHoleCleaner(self, max_size)
        cleaner.run()
        logging.info("Done cleaning holes")

    def set_new_values(self, values):
        self.data.set_new_values(values)

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

                self.data.fill_existing_hole(node_id, start, end, True)

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
        logging.info("Sanitize done")

    def find_valued_areas(self, value):
        logging.info("Finding valued areas equal to %d" % value)
        changes = np.diff(self.data._values == value)
        if value:
            return SparseAreasDict({node_id: self.data.find_valued_areas(node_id, value, changes)
                               for node_id in self.data._graph.blocks
                                }, graph=self.graph)
        else:
            return SparseAreasDict({node_id: self.data.find_valued_areas(node_id, value, changes)
                                   for node_id in self.data._touched_nodes
                                    }, graph=self.graph)
    @classmethod
    def from_intervals(cls, graph, intervals):
        starts, ends = intervals_to_start_and_ends(graph, intervals)
        return cls.from_starts_and_ends(graph, starts, ends)

    @classmethod
    def from_starts_and_ends(cls, graph, starts, ends, dtype=bool):
        pileup = cls(graph)
        for node in starts:
            for i, start in enumerate(starts[node]):
                end = ends[node][i]
                pileup.data.add_value(node, start, end, 1)

        return pileup

    def to_subgraphs(self):
        # Returns a list of areas which each is a subgraph
        collection = SubgraphCollection.from_pileup(self.graph, self)
        return collection

    @classmethod
    def from_valued_areas(cls, graph, valued_areas, touched_nodes = None):
        pileup = cls(graph)

        if touched_nodes is None:
            nodes = graph.blocks
        else:
            nodes = touched_nodes

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

            assert len(starts) == len(ends)
            for start, end in zip(starts, ends):
                pileup.data.add_value(rp, start, end, 1)

        return pileup

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

    def threshold_copy(self, cutoff):
        new_pileup = self.__class__(self.graph)
        new_pileup.data = self.data.threshold_copy(cutoff)
        return new_pileup

    def threshold(self, cutoff):
        self.data.threshold(cutoff)

    @classmethod
    def from_bed_graph(cls, graph, filename):
        raise NotImplementedError()

    @classmethod
    def from_bed_file(cls, graph, filename):
        raise NotImplementedError()

    def to_bed_file(self, filename):
        raise NotImplementedError()

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

    def to_bed_file(self, filename):
        f = open(filename, "w")
        areas = self.find_valued_areas(True)
        for node_id, idxs in areas.items():
            for i in range(len(idxs)//2):
                interval = (node_id, idxs[2*i], idxs[2*i+1])
                f.write("%s\t%s\t%s\t.\t.\t.\n" % interval)
        f.close()
        return filename

    def equals_old_sparse_pileup(self, old_pileup):
        # For tests to pass temporarily
        for node in self.data.nodes():
            indexes, values = self.data.get_sparse_indexes_and_values(node)
            other_indexes = old_pileup.data[node].all_idxs()
            if not np.all(indexes == other_indexes):
                return False

            other_values = old_pileup.data[node].all_values()
            if not np.allclose(values, other_values):
                return False

        return True

    @classmethod
    def create_from_old_sparsepileup(cls, old_pileup):
        # TMP method for tests to pass
        graph = old_pileup.graph
        pileup = DensePileup(graph)

        for node in old_pileup.data:

            values = old_pileup.data[node].all_values()
            indexes = old_pileup.data[node].all_idxs()

            i = 0
            for start, end in zip(indexes[0:-1], indexes[1:]):
                value = values[i]
                pileup.data.set_values(node, start, end, value)
                i += 1

        return pileup

    def set_area_to_value(self, areas, value):
        for node_id in areas.full_areas:
            self.data.set_full_node_value(node_id, value)

        for node_id, internal_intervals in areas.internal_intervals.items():
            assert node_id > 0
            self.data.set_values(node_id, internal_intervals[0], internal_intervals[1], value)

        for node_id, start in areas.starts.items():
            node_size = self.graph.node_size(node_id)
            pileup_end = start
            pileup_start = 0
            if node_id < 0:
                pileup_end = node_size
                pileup_start = node_size - start
                node_id = -node_id

            #print("   Processing start %d for node %d, adding start-end %d-%d" % (start, node_id, pileup_start, pileup_end))
            self.data.set_values(node_id, pileup_start, pileup_end, value)

    def area_is_not_empty(self, areas):
        # Returns true if area is covered by anything anywhere.
        # Areas is Binary cont areas dict
        for node_id in areas.full_areas:
            self.data.add_value_to_full_node(node_id, 1)

        for node_id, internal_intervals in areas.internal_intervals.items():
            assert node_id > 0
            self.data.add_value(node_id, internal_intervals[0], internal_intervals[1], 1)

        for node_id, start in areas.starts.items():
            node_size = self.graph.node_size(node_id)
            pileup_end = start
            pileup_start = 0
            if node_id < 0:
                pileup_end = node_size
                pileup_start = node_size - start
                node_id = -node_id

            #print("   Processing start %d for node %d, adding start-end %d-%d" % (start, node_id, pileup_start, pileup_end))
            self.data.add_value(node_id, pileup_start, pileup_end, 1)

    def add_area(self, areas):
        for node_id in areas.full_areas:
            self.data.add_value_to_full_node(node_id, 1)

        for node_id, internal_intervals in areas.internal_intervals.items():
            assert node_id > 0
            self.data.add_value(node_id, internal_intervals[0], internal_intervals[1], 1)

        for node_id, start in areas.starts.items():
            node_size = self.graph.node_size(node_id)
            pileup_end = start
            pileup_start = 0
            if node_id < 0:
                pileup_end = node_size
                pileup_start = node_size - start
                node_id = -node_id

            #print("   Processing start %d for node %d, adding start-end %d-%d" % (start, node_id, pileup_start, pileup_end))
            self.data.add_value(node_id, pileup_start, pileup_end, 1)

    @classmethod
    def create_from_binary_continous_areas(cls, graph, areas_list):
        pileup = cls(graph, dtype=np.uint8)
        i = 0
        for areas in areas_list:
            if i % 5000 == 0:
                logging.info("Processing read %d" % i)
            i += 1
            pileup.add_area(areas)

        return pileup

    @classmethod
    def from_sparse_files(cls, graph, base_file_name):
        # TODO: Use np.save(...)
        pileup = cls(graph)
        indexes = np.loadtxt(base_file_name + "_indexes.npy", dtype=np.uint32)
        assert np.all(indexes >= 0)
        values = np.loadtxt(base_file_name + "_values.npy")
        touched_nodes = np.loadtxt(
            base_file_name + "_touched_nodes.npy", dtype=np.uint32)

        diffs = np.ediff1d(values, to_begin=[values[0]])
        pileup_vals = pileup.data._values
        pileup_vals[indexes] = diffs
        pileup_vals = np.cumsum(pileup_vals)
        pileup.data._values = pileup_vals
        pileup.data._touched_nodes = set(list(touched_nodes))

        return pileup

    def to_sparse_files(self, file_base_name, truncate_below=0.05):
        vals = self.data._values
        assert np.all(vals >= 0)
        # TODO: Use np.clip(vals, truncate_below, None)
        vals[np.where(vals < truncate_below)] = 0
        # TODO: Use np.nonzero(np.ediff1d(vals, to_begin=vals[0]))
        indexes = np.where(np.ediff1d(vals, to_begin=[vals[0]]) != 0)
        values = vals[indexes]
        np.savetxt(file_base_name + "_touched_nodes.npy",
                   np.array(list(self.data._touched_nodes)))
        np.savetxt(file_base_name + "_indexes.npy", indexes)
        np.savetxt(file_base_name + "_values.npy", values)

        logging.info("Saved p values indexes/values to files")


class DenseControlSample(DensePileup):
    def get_p_dict(self):
        p_value_dict = defaultdict(dict)
        count_dict = defaultdict(int)
        baseEtoTen = np.log(10)
        for node in self.data.nodes():
            for start, end, val in self.data.index_value_pairs(node):
                if val[1] not in p_value_dict[val[0]]:
                    log_p_val = poisson.logsf(val[1], val[0])
                    p_value_dict[val[0]][val[1]] = -log_p_val/baseEtoTen
                p = p_value_dict[val[0]][val[1]]
                count_dict[p] += end-start

        p_value_dict[0.0][0.0] = -1
        count_dict[-1] = 0

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

        new_values = np.apply_along_axis(translation, 1, self.data._values)
        self.data._values = new_values
        self.data.ndim = 1

    def get_scores(self):
        logging.info("Creating p dict")
        self.get_p_dict()
        logging.info("Creating mapping from p-values to q-values")
        self.get_p_to_q_values()
        logging.info("Computing q values")
        self.get_q_values()

    @classmethod
    def from_sparse_control_and_sample(cls, control, sample):
        logging.info("Creating pileup by combining control and sample")
        pileup = cls(sample.graph, ndim=2)
        pileup.data._values[:,0] = control.data._values
        pileup.data._values[:,1] = sample.data._values

        assert np.all(control.data._nodes == sample.data._nodes)

        pileup.data._touched_nodes = sample.data._touched_nodes.union(control.data._touched_nodes)
        return pileup


class QValuesFinder:

    def __init__(self, sample_pileup, control_pileup):
        self.sample = sample_pileup
        self.control = control_pileup

    def get_q_values_pileup(self):

        self.get_p_values()
        self.get_p_to_q_values()
        q_values = self.get_q_values()
        pileup = DensePileup(self.sample.graph)
        pileup.data._values = q_values
        pileup.data._touched_nodes = \
            self.sample.data._touched_nodes.union(self.control.data._touched_nodes)

        return pileup

    def get_p_values(self):
        self.p_values = poisson.logsf(self.sample.data._values,
                                self.control.data._values)
        baseEtoTen = np.log(10)
        self.p_values = -self.p_values / baseEtoTen

        # Only using p-values where sample is not 0
        non_zero_indices = np.nonzero(self.sample.data._values)
        zero_indices = np.where(self.sample.data._values == 0)[0]
        self.p_values[zero_indices] = 0

        sorted_p_values = sorted(self.p_values[non_zero_indices], reverse=True)
        unique, counts = np.unique(sorted_p_values, return_counts=True)
        sorting = np.argsort(-unique)
        self.unique_p_values = unique[sorting]
        self.counts = counts[sorting]

    def get_p_to_q_values(self):
        p_value_counts = self.counts
        p_to_q_values = {}
        rank = 1
        logN = np.log10(sum(p_value_counts))
        pre_q = None
        for i, p_value in enumerate(self.unique_p_values):
            value_count = p_value_counts[i]
            q_value = p_value + (np.log10(rank) - logN)
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
            if x == 0:
                return 0
            return self.p_to_q_values[x]

        trans = np.vectorize(translation, otypes=[np.float])
        new_values = trans(self.p_values)  # np.apply_along_axis(translation, 0, self.p_values)
        assert len(new_values) == len(self.sample.data._values)

        return new_values

