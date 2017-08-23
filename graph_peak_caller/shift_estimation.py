import pyvg as vg
import re
import pyvg.util
import subprocess
from offsetbasedgraph import IntervalCollection


def write_linear_interval_to_bed_file(file_handler, interval):

    strand = "+"
    if interval.direction == -1:
        strand = "-"
    file_handler.write("%s\t%d\t%d\t.\t0\t%s\n" % (interval.start_position.region_path_id,
                                       interval.start_position.offset,
                                       interval.end_position.offset,
                                       strand))


def get_shift_size_on_offset_based_graph(offset_based_graph, interval_file_name):
    linear_graph, trans_to_linear = offset_based_graph.get_arbitrary_linear_graph()
    bed_file = open("tmp2.bed", "w")

    for interval in IntervalCollection.create_generator_from_file(interval_file_name):
        interval.graph = offset_based_graph
        linear_intervals = trans_to_linear.translate_interval(interval).get_single_path_intervals()
        assert len(linear_intervals) == 1
        linear_interval = linear_intervals[0]
        if linear_interval != interval:
            write_linear_interval_to_bed_file(bed_file, linear_interval)


    bed_file.close()
    command = ["macs2", "predictd", "-i", "tmp2.bed", "-g", str(linear_graph.number_of_basepairs()), "-m", "1", "50"]
    output = subprocess.check_output(command, stderr=subprocess.STDOUT)
    output = output.decode("utf-8")

    if "Too few paired peaks" in str(output):
        print(output)
        raise Exception("Too few paired peaks for doing shift estimation")

    tag_size = re.search("tag size = ([0-9]+)", output).groups()[0]
    tag_size = int(tag_size)
    fragment_length = re.search("fragment length is ([0-9]+) bp", output).groups()[0]
    fragment_length = int(fragment_length)
    shift = fragment_length - tag_size
    print("Shift: %d" % shift)
    return fragment_length


if __name__ == "__main__":
    chromosome = "chr4"
    vg_graph = vg.Graph.create_from_file("dm_test_data/x_%s.json" % chromosome, 30000, chromosome)
    #ofbg = vg_graph.get_offset_based_graph()
    print("Intervals:")
    for interval in pyvg.util.vg_mapping_file_to_interval_list(vg_graph, "dm_test_data/reads3_large.json"):
        print("Length: %d" % interval.length())
    print("done")
    #interval_collection = IntervalCollection(intervals)
    #interval_collection.to_file("interval_collection.tmp")
    #get_shift_size_on_offset_based_graph(ofbg, "interval_collection.tmp")

