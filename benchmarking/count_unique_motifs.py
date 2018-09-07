

def main(filename):
    ids = {line.split()[2] for line in open(filename) if not line.startswith("#")}
    return len(ids)


if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
    print(main(filename))
