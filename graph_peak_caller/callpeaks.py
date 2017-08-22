from offsetbasedgraph import IntervalCollection

class CallPeaks(object):
    def __init__(self, graph_file_name, sample_file_name, control_file_name=None, chromosome=None):
        self.graph_file_name = graph_file_name
        self.sample_file_name = sample_file_name
        self.control_file_name = control_file_name if control_file_name is not None else sample_file_name
        self.has_control = control_file_name is not None
        self.linear_genome_size = 0  # Size of linear genome
        self.find_info()
        self.create_graphs()

        self.chromosome = chromosome

    def remove_alignments_not_in_graph(self):
        for alignment_file in [self.sample_file_name, self.control_file_name]:
            if alignment_file is not None:
                interval_collection = IntervalCollection.create_generator_from_file(alignment_file)
                filtered_file = open(alignment_file + "_filtered", "w")
                for interval in self._get_intervals_in_ob_graph(innterval_collection):
                    filtered_file.writelines(["%s\n" % interval.to_file_line()])
                filtered_file.close()

    def _get_intervals_in_ob_graph(self, intervals):
        # Returns only those intervals that exist in vg_graph
        for interval in intervals:
             if interval.region_paths[0] in self.vg_graph.blocks:
                 yield interval

    def find_info(self):
        genome_size = 0
        lines = (line["node"] for line in self.graph_file_name.readlines() if "node" in line)
        sizes = (sum(Node.from_json(json_obj).n_basepairs for json_obj in line) for line in lines)

        self.genome_size = sum(sizes)
        self.n_reads = sum(1 for line in open(self.control_file_name))
        
    def determine_shift(self):
        self.shift = get_shift_size(self.vg_graph, self.sample_file_name, self.chromosome, self.linear_genome_size)
        
    def create_graphs(self):
        self.vg_graph = vg.Graph.create_from_file(self.graph_file_name)
        self.ob_graph = self.vg_graph.get_offset_based_graph()

    def create_control(self):
        bg_track = BackgroundTrack(self.ob_graph, control_file_name, self.d, 1000, 10000, )
        f = open()
        jsons = (j
                 son.loads(line) for line in f.readlines())
        alignments = [vg.Alignment.from_json(json_dict) for json_dict in jsons])
                                   
    def _write_vg_alignments_as_intervals_to_bed_file(self):
        pass

    def write_alignments_to_linear_genome_to_bed_file(self):
        # Write only those alignments that fall on the lniear genome to bed file
        # TODO use method from shift_estimation.py
        pass
