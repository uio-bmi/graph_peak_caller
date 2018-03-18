class Intervals:
    def __init__(self, intervals):
        self._intervals = intervals
        self.n_reads = 0
        self.n_duplicates = 0

    def _count_and_return(self, x):
        self.n_reads += 1
        return x

    def __iter__(self):
        for interval in self._intervals:
            self.n_reads += 1
            yield interval


class UniqueIntervals:
    def __init__(self, intervals):
        self._intervals = intervals
        self.n_reads = 0
        self.n_duplicates = 0

    def __iter__(self):
        interval_hashes = {}
        for interval in self._intervals:
            _hash = interval.hash(ignore_end_pos=True)
            if _hash in interval_hashes:
                self.n_duplicates += 1
                continue

            interval_hashes[_hash] = True
            self.n_reads += 1
            yield interval
