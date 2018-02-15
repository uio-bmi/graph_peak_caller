import numpy as np
from operator import itemgetter


class DiscreteEventSorter(object):
    def __init__(self, index_lists, values_list=None, names=None):
        if names is not None:
            [setattr(self, name.upper(), i) for i, name in enumerate(names)]

        n_types = len(index_lists)
        coded_index_list = [np.array(index_list)*n_types + r
                            for r, index_list in enumerate(index_lists)]

        coded_indices = np.concatenate(coded_index_list)
        sorted_args = np.argsort(coded_indices)
        coded_indices = coded_indices[sorted_args]
        self.indices = coded_indices // n_types
        self.codes = coded_indices % n_types
        if values_list is not None:
            self.values = np.concatenate(values_list)[sorted_args]

    def __iter__(self):
        return zip(self.indices, self.codes, self.values)

    def pileup(self):
        codes = self.codes*2-1
        values = np.cumsum(codes)
        return self.indices, values


class EventSorter(object):

    def np_init(self, index_lists, values_lists):
        all_indices = np.concatenate(index_lists)
        all_values = np.concatenate(values_lists)
        codes = np.empty(all_indices.shape, dtype="uint8")
        start = 0
        for i, index_list in enumerate(index_lists):
            codes[start:start+len(index_list)] = i
            start += len(index_list)
        args = np.lexsort((codes, all_indices))
        self.values = all_values[args]
        self.indices = all_indices[args]
        self.codes = codes[args]
        self.n_codes = len(index_lists)

    def __init__(self, index_lists, values_lists, names=None):

        if names is not None:
            [setattr(self, name.upper(), i) for i, name in enumerate(names)]

        return self.np_init(index_lists, values_lists)

    def __str__(self):
        return "EventSorter(\n" +"\n".join(str(t) for t in self) + ")"

    def __iter__(self):
        # return self.tuples.__iter__()
        return zip(self.indices, self.codes, self.values)


class EventSort(object):
    def __init__(self, index_lists, codes, names=None):
        #print("Running eventnsorter for indexes %s with codes %s" % (index_lists, codes))
        if names is not None:
            [setattr(self, name.upper(), i) for i, name in enumerate(names)]
        self.tuples = []
        for index_list, code in zip(index_lists, codes):
            self.tuples.extend([(idx, code) for idx in index_list])
        self.tuples.sort(key=itemgetter(0, 1))
        values = np.array([t[1] for t in self.tuples])
        indices = np.array([t[0] for t in self.tuples])
        values = np.cumsum(values)
        new_values = np.nonzero(np.diff(values))[0] + 1
        self.values = np.empty(new_values.size+1)
        self.indices = np.empty(new_values.size+1, dtype="int")

        self.indices[1:] = indices[new_values]
        self.indices[0] = indices[0]

        self.values[1:] = values[new_values]
        self.values[0] = values[0]

class DiscreteEventSort(object):
    def __init__(self, index_lists, codes, names=None):
        if names is not None:
            [setattr(self, name.upper(), i) for i, name in enumerate(names)]
        self.tuples = []
        for index_list, code in zip(index_lists, codes):
            self.tuples.extend([(idx, code) for idx in index_list])
        self.tuples.sort(key=itemgetter(0, 1))
        values = np.array([t[1] for t in self.tuples])
        indices = np.array([t[0] for t in self.tuples])
        values = np.cumsum(values)
        new_values = np.nonzero(np.diff(values))[0] + 1
        self.values = np.empty(new_values.size+1)
        self.indices = np.empty(new_values.size+1, dtype="int")

        self.indices[1:] = indices[new_values]
        self.indices[0] = indices[0]
        self.values[1:] = values[new_values]
        self.values[0] = values[0]
