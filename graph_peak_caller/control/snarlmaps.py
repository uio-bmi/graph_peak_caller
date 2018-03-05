import pickle
import logging
import numpy as np

from ..densepileup import DensePileup
from ..sparsediffs import SparseDiffs
from .linearintervals import LinearIntervalCollection


class LinearSnarlMapNP:
    def __init__(self, starts_ends, ):
        self._start_ends = start_ends


class LinearSnarlMap(object):
    def __init__(self, starts, ends, length, graph):
        self._graph = graph
        self._length = length
        self._linear_node_starts, self._linear_node_ends = starts, ends

    @classmethod
    def from_snarl_graph(cls, snarl_graph, graph):
        length = snarl_graph.length()
        starts, ends = snarl_graph.get_distance_dicts()
        return cls(starts, ends, length, graph)

    def get_node_start(self, node_id):
        return self._linear_node_starts[node_id]

    def __str__(self):
        out = "Linear snarl map \n"
        out += " Starts: \n"
        for node, val in self._linear_node_starts.items():
            out += "  %d, %.3f \n" % (node, val)
        out += " Ends: \n"
        for node, val in self._linear_node_ends.items():
            out += "  %d, %.3f \n" % (node, val)

        return out

    def __repr__(self):
        return self.__str__()

    def get_node_end(self, node_id):
        return self._linear_node_ends[node_id]

    def get_scale_and_offset(self, node_id):
        linear_length = self.get_node_end(node_id) \
                        - self.get_node_start(node_id)
        node_length = self._graph.node_size(node_id)
        scale = linear_length/node_length
        offset = self.get_node_start(node_id)
        return scale, offset

    def to_sparse_pileup(self, unmapped_indices_dict, min_value=0):
        all_indices = []
        all_values = []
        i = 0
        node_idxs = self._graph.node_indexes
        min_idx = self._graph.min_node
        for i in range(node_idxs.size-1):
            if i % 100000 == 0:
                logging.info("Processing node %d" % i)
            node_id = i+min_idx
            if node_id not in unmapped_indices_dict:
                all_indices.append(node_idxs[i])
                all_values.append(min_value)
                continue
            unmapped_indices = unmapped_indices_dict[node_id]
            scale, offset = self.get_scale_and_offset(node_id)
            new_idxs = [(idx-offset)//scale+node_idxs[i]
                        for idx in unmapped_indices.indices]
            new_idxs[0] = max(node_idxs[i], new_idxs[0])
            # assert new_idxs[0] == node_idxs[i]
            # assert new_idxs[-1] < node_idxs[i+1]
            all_indices.extend(new_idxs)
            all_values.extend(unmapped_indices.values)
        return SparseDiffs(
            np.array(all_indices, dtype="int"),
            np.diff(np.r_[0, all_values]))

    def to_dense_pileup(self, unmapped_indices_dict):
        pileup = DensePileup(self._graph, dtype=np.uint8)
        i = 0
        for node_id, unmapped_indices in unmapped_indices_dict.items():
            if i % 100000 == 0:
                logging.info("Processing node %d" % i)
            i += 1

            scale, offset = self.get_scale_and_offset(node_id)
            new_idxs = (unmapped_indices.get_index_array()-offset) / scale
            new_idxs = new_idxs.astype("int")
            new_idxs[0] = max(0, new_idxs[0])

            length = self._graph.node_size(node_id)
            indexes = new_idxs
            values = unmapped_indices.get_values_array()
            start_value = values[0]
            # Sanitize indexes
            diffs = np.where(np.diff(indexes) > 0)[0]
            indexes = indexes[diffs+1]
            values = values[diffs+1]

            indexes = np.insert(indexes, 0, 0)
            indexes = np.append(indexes, length)
            values = np.insert(values, 0, start_value)

            j = 0
            for start, end in zip(indexes[:-1], indexes[1:]):
                value = values[j]
                pileup.data.set_values(node_id, start, end, value)
                j += 1
        return pileup

    def map_graph_interval(self, interval):
        start_pos = self.graph_position_to_linear(interval.start_position)
        end_pos = self.graph_position_to_linear(interval.end_position)
        return start_pos, end_pos

    def graph_position_to_linear(self, position):
        node_id = abs(position.region_path_id)
        node_start = self._linear_node_starts[node_id]
        node_end = self._linear_node_ends[node_id]
        node_size = self._graph.node_size(position.region_path_id)
        scale = (node_end-node_start) / node_size
        if position.region_path_id > 0:
            return node_start + scale*position.offset
        else:
            return node_end - scale*position.offset

    def map_interval_collection(self, interval_collection):
        starts = []
        ends = []
        for interval in interval_collection:
            start, end = self.map_graph_interval(interval)
            starts.append(start)
            ends.append(end)
        return LinearIntervalCollection(starts, ends)

    def to_json_files(self, base_name):
        with open(base_name+".length", "w") as f:
            f.write("%s" % self._length)
        with open(base_name+"_starts.pickle", "wb") as f:
            pickle.dump(self._linear_node_starts, f)
        with open(base_name+"_ends.pickle", "wb") as f:
            pickle.dump(self._linear_node_ends, f)

    @classmethod
    def from_json_files(cls, base_name, graph):
        with open(base_name+".length") as f:
            length = int(f.read().strip())
        with open(base_name+"_starts.pickle", "rb") as f:
            starts = pickle.loads(f.read())
        with open(base_name+"_ends.pickle", "rb") as f:
            ends = pickle.loads(f.read())
        return cls(starts, ends, length, graph)

    def to_file(self, file_name):
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def from_file(cls, file_name):
        with open("%s" % file_name, "rb") as f:
            obj = pickle.loads(f.read())
            assert isinstance(obj, cls)
            return obj
