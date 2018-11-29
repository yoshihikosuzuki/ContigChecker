import pickle
import numpy as np
import pandas as pd
from logzero import logger


# Sum up start/end count around peaks (aggregated counts disappear except the peaks)
def aggregate_counts(a, count_threshold, window_size):
    x_spike = np.where(a >= count_threshold)[0]
    x_done = np.zeros(len(a))
    for x in x_spike:
        x_done[x] = a[x]   # ex. [0, 0, ..., 0, 5, 0, ..., 0, 8, 0, ..., 0, 0]
    while sum(x_done) > 0:
        max_pos = np.argmax(x_done)
        s = 0
        for i in range(max_pos - int(window_size / 2), max_pos + int(window_size / 2) + 1):
            s += a[i]
            a[i] = 0
            x_done[i] = 0
        a[max_pos] = s
    return np.where(a >= count_threshold)[0]


def find_spike(a_proper, a_all, count_threshold, window_size):
    x_proper = aggregate_counts(a_proper, count_threshold, window_size)
    x_all = aggregate_counts(a_all, count_threshold, window_size)
    #return sorted(list(set(x_proper) | set(x_all)))
    return (x_proper, x_all)


def annotate_regions():
    # Load files
    contigs = pickle.load(open("contigs.pkl", 'rb'))
    reads = pickle.load(open("reads.pkl", 'rb'))
    mappings = pickle.load(open("mappings.pkl", 'rb'))
    counts = pickle.load(open("counts.pkl", 'rb'))
    depths = pickle.load(open("depths.pkl", 'rb'))

    # Annotate terminal repeats

    # Annotate interspersed repeats

    # Annotate misassemblies
