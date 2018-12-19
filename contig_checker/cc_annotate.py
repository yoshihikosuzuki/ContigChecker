import pickle
import numpy as np
import pandas as pd
from logzero import logger


def find_spikes(count_array,
                count_threshold=3,   # TODO: calculate based on the median depth
                window_size=100):
    pos = np.where(count_array > count_threshold)[0]
    val = count_array[pos]
    spikes = []
    while np.sum(val) > 0:
        # Iteratively find the highest spike and ignore values around the spike
        peak_pos = pos[np.argmax(val)]
        spikes.append(peak_pos)
        val[np.where((peak_pos - int(window_size / 2) < pos)
                     & (pos < peak_pos + int(window_size / 2) + 1))[0]] = 0
    return np.array(spikes)


def find_repeats(x_starts,   # positions of start spikes
                 x_ends,   # end spikes
                 depth):   # depth transition
    # list all possible intervals
    possible_repeat = []   # l, s, t
    min_repeat_len = 1   # >=(length of 1 tracepoint)?
    for start_pos in x_starts:
        for end_pos in x_ends[np.where(x_ends >= start_pos + min_repeat_len - 1)[0]]:
            possible_repeat.append((end_pos - start_pos + 1, start_pos, end_pos))
    possible_repeat.sort()

    # Depth of the intervals in the ascending order of interval length
    # NOTE: we omit previously-detected repeats from the depth calculation of following intervals
    repeats = []
    while len(possible_repeat) > 0:
        l, s, t = possible_repeat[0]
        md = np.mean(depth[s:t+1])
        margin = 10   # TODO: リピートが隣り合う場合にやばい、他のpossible repeatの区間を避ける？
        ld = np.mean(depth[s-margin:s])
        rd = np.mean(depth[t+1:t+1+margin])
        #print(ld, md, rd)
        if 1.5 * ld <= md and 1.5 * rd <= md:   # TODO: optimize parameter
            repeats.append((s, t))
            i = 0
            while i < len(possible_repeat):
                ll, ss, tt = possible_repeat[i]
                if ss == s or tt == t:
                    possible_repeat.pop(i)
                else:
                    i += 1
        else:
            possible_repeat.pop(0)
    return repeats


def annotate_regions():
    ## Load files
    #contigs = pickle.load(open("contigs.pkl", 'rb'))
    #reads = pickle.load(open("reads.pkl", 'rb'))
    #mappings = pickle.load(open("mappings.pkl", 'rb'))
    counts = pickle.load(open("counts.pkl", 'rb'))
    depths = pickle.load(open("depths.pkl", 'rb'))

    # Find start/end spikes
    spikes = (counts[["array_n_start",
                      "array_n_end",
                      "array_n_break_start",
                      "array_n_break_end"]]
              .applymap(lambda df: find_spikes(df)))
    spikes.to_pickle("spikes.pkl")
    #spikes = pickle.load(open("spikes.pkl", 'rb'))

    ## Annotate repeats and misassemblies
    repeats = spikes.apply(lambda df: find_repeats(df["array_n_start"],
                                                   df["array_n_end"],
                                                   depths.loc[df.name, "depth"]), axis=1)
    breaks = spikes.apply(lambda df: find_repeats(df["array_n_break_start"],
                                                  df["array_n_break_end"],
                                                  depths.loc[df.name, "break_depth"]), axis=1)
    annotations = pd.DataFrame({"repeats": repeats, "breaks": breaks})
    annotations.to_pickle("annotations.pkl")
    #annotations = pickle.load(open("annotations.pkl", 'rb'))
