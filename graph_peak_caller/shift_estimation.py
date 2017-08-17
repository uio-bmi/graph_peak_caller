import vg as vg
from vg import vg_mapping_file_to_interval_list

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

    strand = "+"
    if interval.direction == -1:
        strand = "-"
    file_handler.write("%s\t%d\t%d\t%s\n" % (interval.start_position.region_path_id,
                                       interval.start_position.offset,
                                       interval.end_position.offset,
                                       strand))


if __name__ == "__main__":
    chromosome = "chr4"
    vg_graph = vg.Graph.create_from_file("dm_test_data/x_%s.json" % chromosome, 30000, chromosome)
    vg_graph.to_file("%s.vggraph" % chromosome)
    #vg_graph = vg.Graph.from_file("%s.vggraph" % chromosome)

    write_vg_alignemnts_as_linear_intervals_to_bed_file(vg_graph, "dm_test_data/reads3_large.json", chromosome)
