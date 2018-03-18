import logging


class PeakFasta:
    def __init__(self, sequence_retriever):
        self._sequence_retriever = sequence_retriever

    def write_max_path_sequences(self, file_name, max_paths):
        f = open(file_name, "w")
        i = 0
        for max_path in max_paths:
            seq = self._sequence_retriever.get_interval_sequence(max_path)
            f.write(">peak" + str(i) + " " +
                    max_path.to_file_line() + "\n" + seq + "\n")
            i += 1
        f.close()
        logging.info("Wrote max path sequences to fasta file: %s" %
                     (file_name))

    def save_intervals(self, out_fasta_file_name, interval_collection):
        f = open(out_fasta_file_name, "w")
        i = 0
        for max_path in interval_collection.intervals:
            seq = self._sequence_retriever.get_interval_sequence(max_path)
            f.write(">peak" + str(i) + " " +
                    max_path.to_file_line() + "\n" + seq + "\n")
            i += 1
            if i % 100 == 0:
                logging.info("Writing sequence # %d" % i)
        f.close()
