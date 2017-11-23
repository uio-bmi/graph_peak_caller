from offsetbasedgraph import Interval, Block, Translation, GraphWithReversals
import numpy as np
import random
from math import floor
random.seed(1)
from graph_peak_caller.snarls import SimpleSnarl


class SimulatedGraph():
    def __init__(self, graph, translation, n_linear_paths, linear_graph_length, snarls):
        self.graph = graph
        self.translation = translation
        self.n_linear_paths = n_linear_paths
        self.linear_graph_length = linear_graph_length
        self.snarls = snarls

    def translate(self, interval, inverse=False):
        trans = self.translation

        if not inverse:
            interval.graph = trans.graph1
        else:
            interval.graph = self.graph

        return trans.translate_interval(interval, inverse=inverse)\
                    .get_single_path_intervals()[0]


class GraphSimulator():

    def __init__(self, n_paths, n_basepairs_length, n_snps):
        self.n_paths = n_paths
        self.n_snps = n_snps
        self.n_basepairs_length = n_basepairs_length

        self._linear_paths = []

        self.graph = None
        self.translation = None

    def get_simulated_graph(self):
        self._create_graph_with_linear_blocks()
        self._add_snps()
        return SimulatedGraph(self.graph, self.translation,
                              self.n_paths, self.n_basepairs_length,
                              snarls=self.get_snarls_on_snp_graph())

    def _create_graph_with_linear_blocks(self):
        blocks = {i: Block(self.n_basepairs_length)
                  for i in range(100, 100 + self.n_paths)}
        graph = GraphWithReversals(blocks, {})

        # Add dummy blocks at start and end
        start = Block(1)
        end = Block(1)
        for block in graph.blocks:
            graph.adj_list[1].append(block)
            graph.reverse_adj_list[block].append(1)
            graph.adj_list[block].append(2)
            graph.reverse_adj_list[2].append(block)


        graph.blocks[1] = start
        graph.blocks[2] = end

        self.graph = graph
        self.translation = Translation({}, {}, graph)



    def _add_snps(self):
        snp_offsets = 15 * np.array(
                random.sample(range(1, floor((self.n_basepairs_length - 4) / 15)), self.n_snps))
        #print("Snp offsets: " % snp_offsets)
        for snp_offset in snp_offsets:
            #print("Creating snp at %d" % snp_offset)
            paths_to_merge = random.sample(range(0, self.n_paths), 2)
            path1 = paths_to_merge[0] + 100
            path2 = paths_to_merge[1] + 100
            #print("  Paths to merge: %d, %d" % (path1, path2))

            snp_offset = int(snp_offset)

            # Merge before snp
            interval1_before = Interval(snp_offset - 3, snp_offset, [path1], self.graph)
            interval2_before = Interval(snp_offset - 3, snp_offset, [path2], self.graph)

            interval1_before = self.translation.translate(interval1_before)
            interval2_before = self.translation.translate(interval2_before)

            new_graph, trans = self.graph.merge([interval1_before, interval2_before])

            self.graph = new_graph
            self.translation += trans

            # Merge after snp
            interval1_after = Interval(snp_offset + 1, snp_offset + 4, [path1], self.graph)
            interval2_after = Interval(snp_offset + 1, snp_offset + 4, [path2], self.graph)

            interval1_after = self.translation.translate(interval1_after)
            interval2_after = self.translation.translate(interval2_after)

            new_graph, trans = self.graph.merge([interval1_after, interval2_after])

            self.graph = new_graph
            self.translation += trans

        #print(self.graph)
        #print(self.translation)

    def get_snarls_on_snp_graph(self):
        snarls = {}
        max_id = self.graph.max_block_id()


        # Add dummy start and end to graph
        #new_start = max_id + 1
        #new_end = new_start + 1
        #self.graph._add_edge(new_start, 1)
        #self.graph._add_edge(2, new_end)

        #parent_snarl = SimpleSnarl(new_start, new_end, 1, parent=None)
        #snarls[1] = parent_snarl

        current_node = self.graph.get_first_blocks()[0]
        snarl_start = None
        snarl_id = 1000
        while True:
            nexts = self.graph.adj_list[current_node]
            n_nodes_into = len(self.graph.reverse_adj_list[current_node])
            if n_nodes_into == 2:
                assert snarl_start is not None
                snarls[snarl_id] = SimpleSnarl(snarl_start, current_node,
                                               snarl_id) #, parent=parent_snarl)
                #parent_snarl.children.append(snarls[snarl_id])
                snarl_start = None
                snarl_id += 1

            if len(nexts) == 2:
                snarl_start = current_node

            if len(nexts) == 0:
                break

            assert len(nexts) <= 2, "This method only works for snp graphs"
            current_node = nexts[0]

        return snarls


if __name__ == "__main__":
    simulator = GraphSimulator(2, 30, 1)
    simulated_graph = simulator.get_simulated_graph()
    #on_graph = simulated_graph.translate(Interval(2, 20, [100]))
    #print(on_graph)

    print(simulated_graph.snarls)

    #simulator = GraphSimulator(4, 200, 8)





