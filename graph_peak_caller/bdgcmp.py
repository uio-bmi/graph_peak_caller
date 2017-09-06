import shutil
from .pileup import Pileup
import subprocess
import numpy as np


def create_background_pileup_as_max_from_pileups(graph, pileups, background_value, out_file=None, verbose=False):

    """
    first_pileup = pileups.__next__()
    max = first_pileup.get_count_arrays()
    for p in pileups:
        for node_id, count_array in p.get_count_arrays().items():
            max[node_id] = np.maximum(max[node_id], count_array[node_id])

    max_pileup = Pileup(graph)
    max_pileup.set_count_arrays(max)
    if out_file:
        return out_file
    return max_pileup
    """

    if verbose:
        print("Creating bedgraph files")
    pileup_files = [p.to_bed_graph("pileup_%d.tmp" % i)
                    for i, p in enumerate(pileups)]

    max_pileup = "pileup_max.tmp" if out_file is None else out_file

    shutil.copyfile("pileup_0.tmp", max_pileup)

    for pileup in pileup_files:
        if verbose:
            print("Finding max pileup of %s" % pileup)
        max_of_two_pileups(max_pileup, pileup, max_pileup)

    max_pileup_vs_value(max_pileup, background_value, max_pileup)
    if out_file:
        return out_file

    return Pileup.from_bed_graph(graph, max_pileup)


def max_of_two_pileups(pileup1_filename, pileup2_filename, out_filename):
    command = ["macs2", "bdgcmp", "-m", "max", "-t", pileup1_filename, "-c", pileup2_filename, "-o", out_filename]
    output = subprocess.check_output(command)
    #print(output)


def max_pileup_vs_value(pileup_name, value, out_filename):
    command = ["macs2", "bdgopt", "-i" ,pileup_name, "-m", "max", "-p",  str(value),  "-o", out_filename]
    output = subprocess.check_output(command)
    #print(output)


def get_p_value_track(graph, control_file_name, sample_file_name, out_filename):
    command = ["macs2", "bdgcmp", "-t",  sample_file_name, "-c", control_file_name, "-m", "ppois", "-o", out_filename]
    print(" ".join(command))
    output = subprocess.check_output(command)
    return Pileup.from_bed_graph(graph, out_filename)


def get_p_value_track_from_pileups(graph, control_pileup, sample_pileup):
    # Using old macs
    from scipy.stats import poisson

    p_value_pileup = Pileup(graph)
    p_value_count = {}
    for node_id in sample_pileup.get_count_arrays():
        p_value_count[node_id] = -np.log10(1 -
                                    poisson.cdf(sample_pileup.get_count_arrays()[node_id],
                                                control_pileup.get_count_arrays()[node_id])
                                    )

    p_value_pileup.set_count_arrays(p_value_count)
    return p_value_pileup


def create_p_to_q_dict(p_value_counts):
    p_to_q_values = {}
    sorted_p_values = sorted(p_value_counts.keys(), reverse=True)
    rank = 1
    logN = np.log10(sum(p_value_counts.values()))
    pre_q = None
    for p_value in sorted_p_values:
        value_count = p_value_counts[p_value]
        q_value = p_value + (np.log10(rank) - logN)
        if rank == 1:
            q_value = max(0, q_value)
        else:
            q_value = max(0, min(pre_q, q_value))
        p_to_q_values[p_value] = q_value
        pre_q = q_value
        rank += value_count
    return p_to_q_values


def get_q_values_track_from_p_values(p_value_track):
    p_value_counts = p_value_track.count_values()
    p_to_q_values = create_p_to_q_dict(p_value_counts)
    p_value_track.map_values(p_to_q_values)


def scale_down_tracks(ratio, track1, track2):
    if ratio > 1:
        _scale_track(1/ratio, track1)
    else:
        _scale_track(ratio, track2)


def _scale_track(ratio, track):
    command = ["macs2", "bdgopt", "-i", track,
               "-m", "multiply", "-p", str(ratio),
               "-o", track]
    output = subprocess.check_output(command)
    print(output)
