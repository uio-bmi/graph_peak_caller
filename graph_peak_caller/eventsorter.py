import numpy as np
from operator import itemgetter


class DescriteEventSorter(object):
    def __init__(self, index_lists, values_list, names=None):
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
        self.values = np.concatenate(values_list)[sorted_args]

    def __iter__(self):
        return zip(self.indices, self.codes, self.values)


class EventSorter(object):
    def __init__(self, index_lists, values_lists, names=None):
        if names is not None:
            [setattr(self, name.upper(), i) for i, name in enumerate(names)]

        self.tuples = []
        i = 0
        for index_list, values_list in zip(index_lists, values_lists):
            self.tuples.extend([(idx, i, value) for
                                idx, value in zip(index_list, values_list)])
            i += 1
        self.tuples.sort(key=itemgetter(0, 1))

    def __iter__(self):
        return self.tuples.__iter__()


class EventSort(object):
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
        self.indices[0] = 0
        self.values[1:] = values[new_values]
        self.values[0] = values[0]
