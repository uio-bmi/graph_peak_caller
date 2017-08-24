import shutil
from .pileup import Pileup
import subprocess

def create_background_pileup_as_max_from_pileups(graph, pileups, background_value, out_file=None, verbose=False):
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
    output = subprocess.check_output(command)
    #print(output)
    return Pileup.from_bed_graph(graph, out_filename)


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
