import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal, multivariate_normal
import subprocess


def generate_tracks(N):
    tracks = np.zeros((2, N*4))
    tracks[:, :N] = multivariate_normal([1, 1], [[1, 1], [1, 1]], N).transpose()
    tracks[0:, N:2*N] = normal(1, 1, (1, N))
    tracks[1:, 2*N:3*N] = normal(1, 1, (1, N))
    return [tracks[0, :], tracks[1, :]]


def create_replicates(track, n_replicates=2, sigma=0.5):
    replicate_tracks = []
    for i in range(0, n_replicates):
        replicate_tracks.append(track + normal(0, sigma, len(track)))

    return replicate_tracks


def write_track_to_bed(track, file_name):
    lines = []

    for i in range(0, len(track)):
        start_position = 20 * i
        end_position = start_position + 10
        score = track[i]
        lines.append("chr1\t%s\t%s\t.\t0\t.\t%.5f\t-1\t%.5f\t5\n" %
                        (start_position, end_position, score, score))

    file = open(file_name, "w")
    file.writelines(lines)
    file.close()


def run_idr(peak_file_name1, peak_file_name2, out_file_name):
    commmand = ["idr", "--quiet", "--samples", peak_file_name1,
                peak_file_name2, "-o", out_file_name]
    ps = subprocess.Popen(commmand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = ps.communicate()


def forbes(track1, track2):
    observed = sum(track1*track2)
    expected = sum(track1)*sum(track2)/float(len(track1))
    return float(observed)/expected


def idr_bed_to_track_array(bed_file_name, track_length, idr_treshold=540):
    track_array = np.zeros(track_length)
    for line in open(bed_file_name):
        l = line.split()
        start = int(l[1])
        idr_score = int(l[4])
        if idr_score >= idr_treshold:
            #print(idr_score)
            track_position = int(start / 20)
            assert track_position < len(track_array) , \
                "Track position %d from start %d is larger than track length %d" % (track_position, start, len(track_array))
            track_array[track_position] = 1

    return track_array


def noisy_idr(experiments, sigma):
    replicates = [create_replicates(track, sigma=sigma) for track in experiments]

    for name, experiment in zip(["A", "B"], replicates):
        for i, replicate in enumerate(experiment):
            write_track_to_bed(np.exp(replicate), name+str(i) + ".bed")

    for name in ["A", "B"]:
        run_idr(name+"0.bed", name+"1.bed", name + "out.bed")
    final_tracks = {}
    for name in ["A", "B"]:
        final_tracks[name] = idr_bed_to_track_array(name+"out.bed", N*4)

    return forbes(final_tracks["A"], final_tracks["B"]), sum(final_tracks["A"]), sum(final_tracks["B"])


def pipeline(N, sigma):
    experiments = generate_tracks(N)
    return noisy_idr(experiments, sigma)


def test_with_noise(sigma, repeats=100, size=100):
    print("Noise test %.4f" % sigma)
    experiments = generate_tracks(size)
    forbes_values = np.zeros(repeats)
    lengths_A = np.zeros(repeats)
    lengths_B = np.zeros(repeats)
    for i in range(repeats):
        print("repeat %d/%d" % (i, repeats))
        values = noisy_idr(experiments, sigma)
        forbes_values[i] = values[0]
        lengths_A = values[1]
        lengths_B = values[2]

    return np.array([np.mean(forbes_values), np.mean(lengths_A), np.mean(lengths_B)])


N = 20
track_length = N * 4
xs = np.linspace(0.01, 0.4, num=5)
experiment_results = np.array([test_with_noise(x, 400, N) for x in xs])
for i in range(0, len(experiment_results[0])):
    ys = experiment_results[:,i]
    print(ys)
    plt.plot(xs, ys)
    plt.show()
