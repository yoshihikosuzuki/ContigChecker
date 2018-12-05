import pickle
import numpy as np
import pandas as pd
from logzero import logger


# TODO: refactoring
def aggregate_counts(count_array,   # ex. array_n_start
                     count_threshold=3,   # TODO: calculate from median depth
                     window_size=100):
    """
    Find start/end spikes of mapped reads using histogram density estimation with a specific window size.
    """

    x_spike = np.where(count_array >= count_threshold)[0]   # TODO: filtering after density estimation
    x_done = np.zeros(len(count_array), dtype=int)
    for x in x_spike:
        x_done[x] = count_array[x]   # TODO: refactering of this filtering
    while sum(x_done) > 0:
        max_pos = np.argmax(x_done)
        s = 0
        for i in range(max_pos - int(window_size / 2), max_pos + int(window_size / 2) + 1):
            s += count_array[i]
            count_array[i] = 0   # TODO: make it a non-destructive method or copy it first
            x_done[i] = 0
        count_array[max_pos] = s
    return np.where(count_array >= count_threshold)[0]   # array of coordinates of the spikes


def _annotate_regions(counts):


    spike_pos = {}
    spike_pos_list = sorted([(x, data["array_n_start_proper"][x], "start", "proper") for x in x_start_proper] + [(x, array_n_start_all[x], "start", "all") for x in x_start_all] + [(x, data["array_n_end_proper"][x], "end", "proper") for x in x_end_proper] + [(x, array_n_end_all[x], "end", "all") for x in x_end_all])
    for i in range(len(spike_pos_list)):
        spike_pos[i] = spike_pos_list[i]
    spike_pos = pd.DataFrame.from_dict(spike_pos, orient="index", columns=("position", "count", "position_type", "mapping_type"))
    
    spike_pos_all = spike_pos[spike_pos["mapping_type"] == "all"]
    
    #display(spike_pos)
    #display(spike_pos_all)
    
    # list all possible intervals
    possible_repeat = []   # l, s, t
    min_repeat_len = 1   # >=(length of 1 tracepoint)?
    for start_pos in spike_pos_all[spike_pos_all["position_type"] == "start"]["position"]:
        for end_pos in spike_pos_all[spike_pos_all["position_type"] == "end"].pipe(lambda df: df[df["position"] >= start_pos + min_repeat_len - 1])["position"]:
            possible_repeat.append((end_pos - start_pos + 1, start_pos, end_pos))
    possible_repeat.sort()

    # Co-occurence of start and end with a single read
    # 最初からリードを覚えておくのは大変なので、spikeを見つけてからmappingに立ち戻る
    for length, start, end in possible_repeat:
        # NOTE,TODO: 同一リードが同一コンティグの同一領域に複数回マップされることがある(!)ので、setではなくlistそのままで計上（その場合リードのマッピング領域を保持しておく必要がある）？
        # TODO: 上のような例は明らかに異常なので、setにしてリピートアノテーションしなくても良い？？？
        start_reads = set(mappings[mappings["contig_id"] == name].pipe(lambda df: df[start - int(window_size / 2) <= df["contig_start"]]).pipe(lambda df: df[df["contig_start"] <= start + int(window_size / 2)])["read_id"])
        end_reads = set(mappings[mappings["contig_id"] == name].pipe(lambda df: df[end - int(window_size / 2) <= df["contig_end"] - 1]).pipe(lambda df: df[df["contig_end"] - 1 <= end + int(window_size / 2)])["read_id"])
        cooccuring_reads = start_reads & end_reads
        #print(length, start, end, len(start_reads), len(end_reads), len(cooccuring_reads), cooccuring_reads)

    # Depth of the intervals in the ascending order of interval length
    # note that we omit previously-detected repeats from the depth calculation of following intervals
    repeats = []
    while len(possible_repeat) > 0:
        l, s, t = possible_repeat[0]
        md = np.mean(depth_all[s:t+1])
        margin = 10   # TODO: リピートが隣り合う場合にやばい、他のpossible repeatの区間を避ける？
        ld = np.mean(depth_all[s-margin:s])
        rd = np.mean(depth_all[t+1:t+1+margin])
        #print(ld, md, rd)
        if 1.5 * ld <= md and 1.5 * rd <= md:   # TODO: adust param
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
    #print(repeats)
    repeat_starts = set([x[0] for x in repeats])
    repeat_ends = set([x[1] for x in repeats])
    
    # list all possible intervals
    possible_break = []   # l, s, t
    min_break_len = 1   # TODO: or max?
    for end_pos in spike_pos_all[spike_pos_all["position_type"] == "end"]["position"]:
        for start_pos in spike_pos_all[spike_pos_all["position_type"] == "start"].pipe(lambda df: df[df["position"] >= end_pos + min_break_len - 1])["position"]:
            if not (end_pos in repeat_ends and start_pos in repeat_starts):
                possible_break.append((start_pos - end_pos + 1, end_pos, start_pos))
    possible_break.sort()

    # Depth of the intervals in the ascending order of interval length
    # note that we omit previously-detected repeats from the depth calculation of following intervals
    breaks = []
    while len(possible_break) > 0:
        l, s, t = possible_break[0]
        md = np.mean(depth_all[s:t+1])
        margin = 10   # TODO: リピートが隣り合う場合にやばい、他のpossible repeat/breakの区間を避ける？
        ld = np.mean(depth_all[s-margin:s])
        rd = np.mean(depth_all[t+1:t+1+margin])
        #print(ld, md, rd)
        if 1.5 * md <= ld and 1.5 * md <= rd:   # TODO: adust param
            breaks.append((s, t))
            i = 0
            while i < len(possible_break):
                ll, ss, tt = possible_break[i]
                if ss == s or tt == t:
                    possible_break.pop(i)
                else:
                    i += 1
        else:
            possible_break.pop(0)
    #print(breaks)
    
    #annot = determine_annotation(spike_pos, depth_all, depth_proper)
    
    return (repeats, breaks)


def annotate_regions():
    ## Load files
    contigs = pickle.load(open("contigs.pkl", 'rb'))
    reads = pickle.load(open("reads.pkl", 'rb'))
    mappings = pickle.load(open("mappings.pkl", 'rb'))
    counts = pickle.load(open("counts.pkl", 'rb'))
    depths = pickle.load(open("depths.pkl", 'rb'))

    ## Annotate terminal/interspersed repeats
    # Start/end spikes
    spike_proper = (counts.loc[pd.IndexSlice[:, "proper"],   # extract only proper
                               ["array_n_start",
                                "array_n_end",
                                "array_n_break_start",
                                "array_n_break_end"]]
                    .applymap(lambda df: aggregate_counts(df))
                    .reset_index("terminal_type").drop("terminal_type", axis=1))   # remove second index
    
    spike_all = (counts[["array_n_start",
                         "array_n_end",
                         "array_n_break_start",
                         "array_n_break_end"]]
                 .groupby("contig_id")
                 .apply(sum)   # proper + clipped
                 .applymap(lambda df: aggregate_counts(df)))   # TODO: check abort

    ## Annotate misassemblies
