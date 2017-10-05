from .extender import AreasBuilder, Areas
import offsetbasedgraph as obg


class Cleaner(object):
    def __init__(self, pileup, threshold):
        self.graph = pileup.graph
        self.areas = self.get_areas(pileup)
        self.areas_builder = AreasBuilder(self.graph)
        self.get_starts_and_ends_dict(self.areas)
        self.intervals = []
        self.threshold = threshold

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
                self.starts_dict[-node] = int(self.graph.node_size(node)-startends[-2])
                self.ends_dict[node] = int(self.graph.node_size(node)-startends[-2])

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
        return [node_list + [next_node] for next_node in self.cur_adj_list[last_node]
                if next_node in self.starts_dict]

    def run(self):
        self.other_adj_list = self.graph.reverse_adj_list
        for adj_list in [self.graph.adj_list, self.graph.reverse_adj_list]:
            print("Running dir")
            self.cur_adj_list = adj_list
            self.directed_run(adj_list)
            self.other_adj_list = self.cur_adj_list
        self.finalize()
        return self.areas

    def directed_run(self, adj_list):
        node_lists = self.get_init_nodes()
        assert all(node_list[0] in self.ends_dict for node_list in node_lists)
        print(node_lists)
        while node_lists:
            new_list = []
            for node_list in node_lists:
                extensions = self.extend_node_list(node_list)
                should_extend = self.handle_node_list(node_list, extensions)
                if not should_extend:
                    continue
                new_list.extend(self.extend_node_list(node_list))
            node_lists = new_list

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

    def save(self, node_list):
        print("Saving", node_list)
        areas = {node_id: [0, self.starts_dict[node_id]]
                 for node_id in node_list[1:]}
        areas.update({-node_list[0]: [0, self.starts_dict[-node_list[0]]]})
        print(areas)
        self.areas_builder.update(areas)

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
                areas[-node_id].append(self.graph.node_size(node_id)-startend[1])
                areas[-node_id].append(self.graph.node_size(node_id))
        print("Finalized", areas)
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

    def save_old(self, node_list):
        print("Saving: ", node_list)
        start = self.graph.node_size(node_list[0])-self.ends_dict[node_list[0]]
        end = self.starts_dict[node_list[-1]]
        interval = obg.DirectedInterval(
            start, end, node_list, graph=self.graph)
        print(interval)
        self.intervals.append(interval)

    def handle_node_list(self, node_list, extensions):
        assert node_list[0] in self.ends_dict
        print("Handling: ", node_list)
        if node_list[-1] in node_list[1:-1]:  # Loop
            print("Loop")
            self.save(node_list)
            return False
        if extensions:
            print("Extending")
            return True
        length = self.get_length(node_list)
        print("Checking length:", length)
        if length >= self.threshold:
            print("Saving")
            self.save(node_list)

        return False


class HolesCleaner(Cleaner):
    def handle_node_list(self, node_list, extensions):
        assert node_list[0] in self.ends_dict
        if node_list[-1] in node_list[1:-1]:
            print("##### LOOP")
            return False
        length = self.get_length(node_list)
        if length > self.threshold:
            return False
        if len(extensions) == len(self.cur_adj_list[node_list[-1]]):
            return True
        print(node_list[-1], extensions, self.cur_adj_list[node_list[-1]])
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
        return end-start <= self.threshold
