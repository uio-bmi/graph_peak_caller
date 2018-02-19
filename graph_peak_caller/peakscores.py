from collections import deque, defaultdict
import logging
import numpy as np
from .peakcollection import Peak
from .sparsepileupv2 import RpScore
from .sparsepileup import ValuedIndexes


class ScoredPeak(object):
    def __init__(self, peak, scores):
        self._peak = peak
        self._scores = scores
        self._graph = peak.graph

    def __eq__(self, other):
        if self._peak != other._peak:
            return False

        for node, scores in self._scores.items():
            if isinstance(scores, ValuedIndexes):
                scores = RpScore.from_valued_indexes(scores)
            if isinstance(other._scores[node], ValuedIndexes):
                scores = RpScore.from_valued_indexes(other._scores[node])

            if scores != scores:
                return False

        return True

    @staticmethod
    def _from_valued_indexes(valued_indexes, start, end):
        return valued_indexes.get_subset(start, end)

    @classmethod
    def from_peak_and_numpy_pileup(cls, peak, pileup):
        scores = {}
        for node_id in peak.full_areas:
            length = pileup.graph.node_size(node_id)
            node_scores = pileup.data.score(node_id, 0, length)
            scores[node_id] = node_scores
            scores[-node_id] = node_scores

        for node_id, start in peak.starts.items():
            length = pileup.graph.node_size(node_id)
            if node_id < 0:
                node_scores = pileup.data.score(-node_id, length-start, length)
            else:
                node_scores = pileup.data.score(node_id, 0, start)
            scores[node_id] = node_scores

        for node_id, startend in peak.internal_intervals.items():
            node_scores = pileup.data.score(node_id, startend[0], startend[1])

            scores[node_id] = node_scores

        return cls(peak, scores)

    @classmethod
    def from_peak_and_pileup(cls, peak, pileup):
        return cls.from_peak_and_numpy_pileup(peak, pileup)

    def __str__(self):
        return "\n".join("%s: (max: %s, sum: %s)" % (node_id, vi[0], vi[1]) for
                         node_id, vi in self._scores.items())

    __repr__ = __str__

    def _get_adj_list(self):
        starts = self._peak.starts.keys()
        ends = [-node_id for node_id in self._peak.starts.keys()]
        fulls = list(self._peak.full_areas.keys())
        adj_list = {node_id:
                    [next_node for next_node in self._graph.adj_list[node_id]
                     if abs(next_node) in fulls or next_node in starts]
                    for node_id in ends+fulls+[-n for n in fulls]}
        d_adj_list = defaultdict(list)
        d_adj_list.update(adj_list)
        return d_adj_list

    @staticmethod
    def __clean_sums(sums):
        values = [s for s in sums.values() if s != np.inf]
        if len(values) > 0:
            max_finite_value = max(values)
        else:
            max_finite_value = 100

        for node_id, sum_value in sums.items():
            if sum_value == np.inf:
                sums[node_id] = max_finite_value+1

    def clean_path(self, path, sums):
        if sums[-(path.region_paths[0])]:
            return path
        for i, region_path in enumerate(path.region_paths[1:]):
            if sums[region_path]:
                offset = np.nonzero(self._scores[region_path])[0][0]
                p = Peak(
                    offset, path.end_position.offset,
                    path.region_paths[i+1:], graph=self._graph)
                p.set_score(path.score)
                return p

    def get_max_path(self):
        logging.info("Processing peak")
        sums = {node_id: float(scores.sum()) for node_id, scores
                in self._scores.items()}
        logging.info(" Sums fetched")
        # Handle peaks that are on one node
        if len(self._peak.internal_intervals) > 0:
            node, start_end = list(self._peak.internal_intervals.items())[0]

            interval = Peak(start_end[0], start_end[1],
                            [node], graph=self._graph)
            score = sums[node]
            interval.set_score(
                self._scores[node].max_value()) # score / interval.length())
            return interval

        logging.info("Clean sums")
        self.__clean_sums(sums)
        for key in list(sums.keys()):
            if key in self._peak.full_areas:
                sums[-key] = sums[key]

        adj_list = self._get_adj_list()
        start_positions = self._peak.get_start_positions()
        ends = [pos.region_path_id for pos in start_positions]
        start_values = [sums[-node_id] if (-node_id in self._peak.starts)
                        else sums[abs(node_id)]
                        for node_id in ends]

        memo = defaultdict(lambda: -1)
        stack = deque(zip([[e] for e in ends], start_values))
        if not any(sums.values()):
            logging.warning(
                "Peak has only values 0: %s",
                self._peak)
            return None
        assert stack, str(self._peak)
        global_max = -1
        global_max_path = None
        n_iterations = 0
        while stack:
            node_ids, value = stack.popleft()

            if memo[node_ids[-1]] >= value:
                continue
            if value >= global_max:
                global_max_path = node_ids
                global_max = value

            memo[node_ids[-1]] = value
            nexts = adj_list[node_ids[-1]]
            new_items = [
                (node_ids+[next_node], value+sums[next_node])
                for next_node in nexts
                if next_node not in node_ids[1:] and next_node in sums]
            stack.extend(new_items)

            if n_iterations % 20 == 0:
                logging.info("  %d iterations" % n_iterations)
            n_iterations += 1

        start_node = global_max_path[0]
        start_pos = [pos for pos in start_positions if
                     pos.region_path_id == start_node][0]
        start_offset = int(start_pos.offset)
        # start_offset = int(self._graph.node_size(start_node) -
        # self._peak.starts[-start_node])
        end_node = global_max_path[-1]
        end_offset = self._peak.starts[end_node] if end_node in self._peak.starts else self._peak.graph.node_size(end_node)
        if -global_max_path[0] in self._scores:
            max_score_in_peak = self._scores[-global_max_path[0]].max_value()
        else:
            max_score_in_peak = self._scores[global_max_path[0]].max_value()

        for node in global_max_path[1:-1]:
            max_in_node = self._scores[abs(node)].max_value()
            max_score_in_peak = max(max_score_in_peak, max_in_node)

        if len(global_max_path) > 1:
            max_score_in_peak = max(
                max_score_in_peak,
                self._scores[global_max_path[-1]].max_value())

        max_path_peak = Peak(
            int(start_offset), int(end_offset),
            global_max_path, graph=self._graph)

        score = max_score_in_peak  # global_max / max_path_peak.length()
        # score = global_max / max_path_peak.length()
        max_path_peak.set_score(score)
        return max_path_peak
