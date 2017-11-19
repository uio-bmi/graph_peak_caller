class CollectionIO(object):
    _obj_type = None

    def __init__(self, collection):
        self._collection = collection

    @classmethod
    def from_file(cls, file_name, graph):
        with open(file_name) as f:
            collection = [cls._obj_type.from_file_line(line, graph)
                          for line in f.readlines()]
        return cls(collection)

    def to_file(self, file_name):
        with open(file_name, "w") as f:
            f.writelines(item.to_file_line() for item in self._collection)
