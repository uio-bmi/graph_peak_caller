
class GraphNotFoundException(Exception):
    pass

class InvalidVgJsonGamFile(Exception):
    def __init__(self, msg=None):
        if msg is None:
            msg = """Vg alignment file is invalid and could not be parsed.
            When input file with .json-ending is used, Graph Peak Caller assumes a vg
            alignment file in json formatet, created by converting a vg .gam file
            by using vg view -aj.
            """
            super(InvalidVgJsonGamFile, self).__init__(msg)

class InvalidPileupInterval(Exception):
    pass