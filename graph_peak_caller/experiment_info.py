import pickle


class ExperimentInfo(object):
    def __init__(self, genome_size, fragment_length, read_length):
        self.genome_size = genome_size
        self.n_sample_reads = 0  # Counters will be modified by Callpeaks
        self.n_control_reads = 0
        self.fragment_length = fragment_length
        self.read_length = read_length

    def to_file(self, file_name):
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def from_file(cls, file_name):
        with open("%s" % file_name, "rb") as f:
            data = pickle.loads(f.read())
            return data
