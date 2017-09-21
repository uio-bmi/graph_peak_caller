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
    n_not_translation = 0

    if not isinstance(interval_file_name, str):
        # We have generator of intervals
        file_name = interval_file_name.to_file("tmp_intervals_for_shifting")
        intervals = IntervalCollection.from_file(file_name)
        interval_file_name = IntervalCollection.from_file(file_name)  # Get back generator
    else:
        intervals = IntervalCollection.from_file(interval_file_name)

    for interval in intervals:
        interval.graph = offset_based_graph
        if trans_to_linear.interval_has_translation(interval):
            linear_intervals = trans_to_linear.translate_interval(interval).get_single_path_intervals()
            assert len(linear_intervals) == 1
            linear_interval = linear_intervals[0]
            linear_interval.direction = interval.direction
            write_linear_interval_to_bed_file(bed_file, linear_interval)
        else:
            n_not_translation += 1

    print("%d intervals did not have translation to linear graph" % n_not_translation)

    bed_file.close()
    command = ["macs2", "predictd", "-i", "tmp2.bed", "-g", str(linear_graph.number_of_basepairs()), "-m", "5", "50"]
    command_str = ' '.join(command)
    print("Macs shift command: " + command_str)
    output = subprocess.check_output(command, stderr=subprocess.STDOUT)
    output = output.decode("utf-8")
    #print(output)
    if "Too few paired peaks" in str(output):
        #print(output)
        raise RuntimeError("Too few paired peaks for doing shift estimation")

    tag_size = re.search("tag size = ([0-9]+)", output).groups()[0]
    tag_size = int(tag_size)
    fragment_length = re.search("fragment length is ([0-9]+) bp", output).groups()[0]
    fragment_length = int(fragment_length)
    shift = fragment_length - tag_size
    print("Shift: %d" % shift)
    return fragment_length, tag_size


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

