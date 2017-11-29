from .extender import AreasBuilder, Areas
import logging
LOOP_VALUE = 3000000000


class Cleaner(object):
    def __init__(self, pileup, threshold):
        self.graph = pileup.graph
        self.areas = self.get_areas(pileup)
        self.areas_builder = AreasBuilder(self.graph)
        self.get_starts_and_ends_dict(self.areas)
        self.intervals = []
        self.threshold = threshold
        self.ignored_nodes = set([])

    def get_starts_and_ends_dict(self, areas):
        self.starts_dict = {}
        self.ends_dict = {}
        for node, startends in areas.items():
            if not startends:
                continue
            node = int(node)
            if startends[0] == 0:
                self.starts_dict[node] = int(startends[1])
                self.ends_dict[-node] = int(startends[1])
            if startends[-1] == self.graph.node_size(node):
                self.starts_dict[-node] = int(
                    self.graph.node_size(node) - startends[-2])
                self.ends_dict[node] = int(
                    self.graph.node_size(node) - startends[-2])

        logging.debug("N starts: %s", len(self.starts_dict))
        logging.debug("N ends: %s", len(self.ends_dict))

    def is_end_included(self, node_id):
        pass

    def is_start_included(self, node_id):
        pass

    def handle_infinite(self, node_list):
        pass

    def extend_node_list(self, node_list):
        last_node = node_list[-1]
        if len(node_list) > 1 and not self._is_region_path_covered(last_node):
            return []
        extended = [node_list + [next_node] for next_node
                    in self.cur_adj_list[last_node]
                    if next_node in self.starts_dict]

        for added_node in extended:
            self.ignored_nodes.discard(added_node[-1])
            self.ignored_nodes.discard(-added_node[-1])
            logging.info("Discarding node %d" % added_node[-1])

        return extended

    def run(self):
        self.other_adj_list = self.graph.reverse_adj_list
        for adj_list, name in zip([self.graph.adj_list, self.graph.reverse_adj_list], ("F", "B")):
            self.cur_dir = name
            self.cur_adj_list = adj_list
            self.directed_run(adj_list)
            self.other_adj_list = self.cur_adj_list
        self.finalize()
        return self.areas

    def directed_run(self, adj_list):
        self._cur_memo = {}
        self._cur_remain_memo = {}
        node_lists = self.get_init_nodes()
        assert all(node_list[0] in self.ends_dict for node_list in node_lists)
        while node_lists:
            logging.debug("N lists: %s", len(node_lists))
            node_list = node_lists.pop()
            extensions = self.extend_node_list(node_list)
            should_extend = self.handle_node_list(node_list, extensions)
            if not should_extend:
                continue
            node_lists.extend(self.extend_node_list(node_list))

    def handle_node_list(self, node_list, extensions):
        raise NotImplementedError

    def get_length(self, node_list):
        length = sum(self.graph.node_size(node_id)
                     for node_id in node_list[1:-1])
        length += self.ends_dict[node_list[0]]
        if len(node_list) > 1:
            length += self.starts_dict[node_list[-1]]
        return length

    def _is_region_path_covered(self, node_id):
        return node_id in self.starts_dict and self.starts_dict[node_id] == self.graph.node_size(node_id)

    def get_init_nodes(self):
        return [[node] for node in self.ends_dict.keys() if
                self.is_init_node(node)]

    def save(self, node_list, memo_value=None):
        areas = {node_id: [0, self.starts_dict[node_id]]
                 for node_id in node_list[1:]}
        areas.update({-node_list[0]: [0, self.starts_dict[-node_list[0]]]})
        self.areas_builder.update(areas)
        if memo_value is not None:
            for node_id in node_list[1:]:
                self._cur_remain_memo[node_id] = memo_value

    def finalize(self):
        areas = {}
        for node_id, startends in self.areas.items():
            new_start_ends = []
            for i in range(len(startends) // 2):
                start = int(startends[i*2])
                end = int(startends[i*2+1])
                if self._check_internal_interval(node_id, start, end):
                    new_start_ends.extend([start, end])
            if new_start_ends:
                areas[node_id] = new_start_ends

        for node_id, startend in self.areas_builder.areas.items():
            if abs(node_id) not in areas:
                areas[abs(node_id)] = []
            if node_id > 0:
                areas[node_id].insert(0, startend[1])
                areas[node_id].insert(0, 0)
            else:
                areas[-node_id].append(
                    self.graph.node_size(node_id)-startend[1])
                areas[-node_id].append(self.graph.node_size(node_id))

        for ignored_node in self.ignored_nodes:
            logging.warning("Ignored node: %d" % ignored_node)
            # areas[abs(ignored_node)] = [0, self.graph.node_size(ignored_node)]
        self.areas = Areas(self.graph, areas)


class PeaksCleaner(Cleaner):
    def is_init_node(self, node):
        if self.ends_dict[node] < self.graph.node_size(node):
            return True
        return not bool([prev_node for prev_node in self.other_adj_list[-node]
                         if -prev_node in self.ends_dict])

    def get_areas(self, pileup):
        return pileup.find_valued_areas(True)

    def _check_internal_interval(self, node_id, start, end):
        if start == 0 or end == self.graph.node_size(node_id):
            return False
        return end-start >= self.threshold

    def handle_node_list(self, node_list, extensions):
        assert node_list[0] in self.ends_dict
        if node_list[-1] in node_list[1:-1]:  # Loop
            self.save(node_list, memo_value=LOOP_VALUE)
            return False

        length = self.get_length(node_list)
        last_node = node_list[-1]
        if last_node in self._cur_memo and length <= self._cur_memo[last_node]:
            if last_node not in self._cur_remain_memo:
                return False
            my_remain = self._cur_remain_memo[last_node] - self._cur_memo[last_node] + length
            if my_remain >= self.threshold:
                self.save(node_list, memo_value=length)
            return False
        self._cur_memo[last_node] = length
        if extensions:
            return True
        length = self.get_length(node_list)
        if length >= self.threshold:
            self.save(node_list, memo_value=length)
        return False

    def get_init_nodes(self):
        init_nodes = []
        for node in self.ends_dict.keys():
            if self.is_init_node(node):
                init_nodes.append([node])
            else:
                if node in self.ends_dict:
                    self.ignored_nodes.add(node)
        return init_nodes


class HolesCleaner(Cleaner):
    def handle_node_list(self, node_list, extensions):
        assert node_list[0] in self.ends_dict
        if node_list[-1] in node_list[1:-1]:
            return False
        length = self.get_length(node_list)
        if length > self.threshold:
            return False

        if (len(self.cur_adj_list[node_list[-1]])) == 0:
            if self._is_region_path_covered(node_list[-1]):
                return False
            if node_list[-1] in self.starts_dict and len(node_list) > 1:
                self.save(node_list)
                return True
            return False

        if len(extensions) and (len(extensions) == len(self.cur_adj_list[node_list[-1]])):
            return True

        self.save(node_list)
        return True

    def get_areas(self, pileup):
        return pileup.find_valued_areas(False)

    def is_init_node(self, node):
        if self.ends_dict[node] < self.graph.node_size(node):
            return True
        return not all([-prev_node in self.ends_dict
                        for prev_node in self.other_adj_list[-node]])

    def _check_internal_interval(self, node_id, start, end):
        if start == 0 or end == self.graph.node_size(node_id):
            return False

        keep = end-start <= self.threshold
        return keep
