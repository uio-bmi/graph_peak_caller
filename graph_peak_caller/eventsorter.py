import numpy as np


class EventSorter(object):
    def __init__(self, index_lists, values_list, names=None):
        if names is not None:
            [setattr(self, name.upper(), i) for i, name in enumerate(names)]

        n_types = len(index_lists)
        coded_index_list = [np.array(index_list)*n_types + r
                            for r, index_list in enumerate(index_lists)]
        coded_indices = np.concatenate(coded_index_list)
        sorted_args = np.argsort(coded_indices)
        coded_indices = coded_indices[sorted_args]
        self.indices = coded_indices//n_types
        self.codes = coded_indices % n_types
        self.values = np.concatenate(values_list)[sorted_args]

    def __iter__(self):
        return zip(self.indices, self.codes, self.values)
