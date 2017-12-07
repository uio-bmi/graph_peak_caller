from graph_peak_caller.eventsorter import EventSorter


def test_event_sorter():
    indices = [[2, 10, 4], [2, 5, 9, ]]
    values = [[1, 2, 3], [5, 6, 7]]
    sorter = EventSorter(indices, values)
    assert list(sorter) == [(2, 0, 1),
                            (2, 1, 5),
                            (4, 0, 3),
                            (5, 1, 6),
                            (9, 1, 7),
                            (10, 0, 2)]
