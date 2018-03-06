from collections import defaultdict


class FimoEntry:
    def __init__(self, peak_id, start, stop, direction, score):
        self.peak_id = peak_id
        self._start = start
        self._end = stop
        self._direction = direction
        self._score = score

    @classmethod
    def from_file_line(cls, line):
        parts = line.split()
        peak_id = parts[2]
        start = int(parts[3])
        stop = int(parts[4])
        direction = parts[5]
        score = float(parts[6])
        return cls(peak_id, start, stop, direction, score)

    def __repr__(self):
        return "FE(%s: %s-%s [%s])" % (
            self.peak_id, self._start, self._end, self._direction)


class FimoFile:
    def __init__(self, entry_dict):
        self._entry_dict = entry_dict

    def get_entries(self, peak_id):
        return self._entry_dict[peak_id]

    @classmethod
    def from_file(cls, filename):
        entries = defaultdict(list)
        for line in open(filename):
            if line.startswith("#"):
                continue
            fimo_entry = FimoEntry.from_file_line(line)
            entries[fimo_entry.peak_id].append(fimo_entry)
        return cls(entries)
