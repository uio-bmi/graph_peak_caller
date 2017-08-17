import vg as vg
from vg import vg_mapping_file_to_interval_list

def estimate_shift(vg_graph, alignments_json_file, limit_to_chromosome):

    print("Creating translation")
    translation = vg_graph.get_translation(limit_to_chromosome)

    print("Mapping to intervals")


    linear_intervals = translate_intervals_to_linear_genome(mappings, translation)
    print("Linear intervals: " , len(linear_intervals))
    intervals_to_bed_file(linear_intervals, "tmp.bed")


def write_vg_alignemnts_as_linear_intervals_to_bed_file(vg_graph, alignments_json_file, limit_to_chromosome):

    bed_file = open("tmp.bed", "w")

    print("Creating translation")
    translation = vg_graph.get_translation(limit_to_chromosome)

    print("Mapping to intervals")
    i = 0
    for alignment_interval in vg_mapping_file_to_interval_list(vg_graph, alignments_json_file):

        #print("Alignment %d" % i)
        # Translate to linear genome
        linear_interval = translate_alignment_to_linear_genome(alignment_interval, translation)
        if linear_interval:
            write_linear_interval_to_bed_file(bed_file, linear_interval)

        i += 1
    bed_file.close()


def translate_alignment_to_linear_genome(alignment_interval, translation):
    if all([rp in translation._b_to_a for rp in alignment_interval.region_paths]):
        translated_intervals = translation.translate_interval(alignment_interval, True).get_single_path_intervals()
        #print(translated_intervals)
        assert len(translated_intervals) <= 1

        if len(translated_intervals) > 0:
            return translated_intervals[0]

    return False

def write_linear_interval_to_bed_file(file_handler, interval):
    file_handler.write("%s\t%d\t%d\n" % (interval.start_position.region_path_id,
                                       interval.start_position.offset,
                                       interval.end_position.offset))

def translate_intervals_to_linear_genome(intervals, translation):
    # Assuming translation is a translation from main path to interval's graph
    filtered_intervals = []
    for interval in intervals:

        if all([rp in translation._b_to_a for rp in interval.region_paths]):

            translated_intervals = translation.translate_interval(interval, True).get_single_path_intervals()
            #print(translated_intervals)
            assert len(translated_intervals) <= 1

            if len(translated_intervals) > 0:
                filtered_intervals.append(translated_intervals[0])

    return filtered_intervals

def intervals_to_bed_file(intervals, file_name):
    print("Creating bed file")
    lines = []
    outfile = open(file_name, "w")
    i = 0
    for interval in intervals:
        #print(interval)
        print("Interval %d / %d" % (i, len(intervals)))
        i += 1
        assert interval.start_position.region_path_id == interval.end_position.region_path_id
        lines.append("%s\t%d\t%d\n" % (interval.start_position.region_path_id,
                                       interval.start_position.offset,
                                       interval.end_position.offset))

    outfile.writelines(lines)
    outfile.close()


if __name__ == "__main__":
    chromosome = "chr4"
    vg_graph = vg.Graph.create_from_file("dm_test_data/x_%s.json" % chromosome, 30000, chromosome)
    vg_graph.to_file("%s.vggraph" % chromosome)
    #vg_graph = vg.Graph.from_file("%s.vggraph" % chromosome)

    write_vg_alignemnts_as_linear_intervals_to_bed_file(vg_graph, "dm_test_data/reads3_large.json", chromosome)
    #estimate_shift(vg_graph, "dm_test_data/reads3_large.json", chromosome)
