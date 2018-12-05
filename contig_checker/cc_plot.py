import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go

plt.style.use('ggplot')
py.init_notebook_mode()

# Visualization tools

# TODO: scale to Mbp size (use window?)
# TODO: vis. of mapped reads

def plot_depth(counts, name, count_threshold=3, window_size=100, width=1000, height=600, out_html=None):
    if not isinstance(name, int):
        name = counts[counts["contig_header"] == name].index[0]
        
    data = counts.loc[name]
    contig_header = data["contig_header"]
    
    ## Depth plot
    depth_all, depth_proper = calc_depth(data)
    
    trace1 = go.Scatter(x=[0, data["contig_length"] - 1],
                       y=[data["n_contains"], data["n_contains"]],
                       mode="lines",
                       line=dict(color="grey"),
                       name="contains")
    
    trace2 = go.Scatter(x=np.arange(data["contig_length"]),
                            y=depth_proper,
                            mode='lines',
                        line=dict(color="dodgerblue"),
                            name="depth (proper)")
    
    trace3 = go.Scatter(x=np.arange(data["contig_length"]),
                            y=depth_all,
                            mode='lines',
                        line=dict(color="limegreen"),
                            name="depth (all)")
    
    trace4 = go.Scatter(x=[0, 0, data["contig_length"] - 1, data["contig_length"] - 1],
                       y=[data["n_dovetail_start_proper"], data["n_dovetail_start_proper"] + data["n_dovetail_start_clipped"], data["n_dovetail_end_proper"], data["n_dovetail_end_proper"] + data["n_dovetail_end_clipped"]],
                       mode="markers",
                       marker=dict(color=["blue", "green", "blue", "green"], symbol=24, size=14),
                       name="dovetail")
    
    ## Start/end count plot
    
    # TODO: 最初に1塩基単位でcount_threshold以上のcountを見つけるのではなく、windowでcountする
    
    # TODO: functionalize. use multiindex of columns ("depth array" -> "depth", "break depth" -> "start", "end")
    
    #count_threshold = 3   # only on the points with > count_threshold start [end] count, plot start [end] data
    array_n_start_all = data["array_n_start_proper"] + data["array_n_start_clipped"]   # TODO: proper>0でall=0のことがある。aggregationの後に再度合計？
    x_start_proper, x_start_all = find_spike(data["array_n_start_proper"], array_n_start_all, count_threshold, window_size)
    x_start = sorted(list(set(x_start_proper) | set(x_start_all)))
    array_n_end_all = data["array_n_end_proper"] + data["array_n_end_clipped"]
    x_end_proper, x_end_all = find_spike(data["array_n_end_proper"], array_n_end_all, count_threshold, window_size)
    x_end = sorted(list(set(x_end_proper) | set(x_end_all)))    
    
    trace5 = go.Scatter(x=x_start,
                    y=data["array_n_start_proper"][x_start],
                        mode="markers",
                        marker=dict(color="blue", symbol=8, size=14),
                   name="start (proper)")
    
    trace6 = go.Scatter(x=x_start,
                    y=array_n_start_all[x_start],
                        mode="markers",
                        marker=dict(color="green", symbol=8, size=14),
                   name="start (all)")
    
    trace7 = go.Scatter(x=x_end,
                    y=data["array_n_end_proper"][x_end],
                        mode="markers",
                        marker=dict(color="blue", symbol=7, size=14),
                   name="end (proper)")
    
    trace8 = go.Scatter(x=x_end,
                    y=array_n_end_all[x_end],
                        mode="markers",
                        marker=dict(color="green", symbol=7, size=14),
                   name="end (all)")
    
    ## Break depth plot
    depth_break_all, depth_break_proper = calc_break_depth(data)
    
    trace9 = go.Scatter(x=np.arange(data["contig_length"]),
                            y=-depth_break_proper,
                            mode='lines',
                        line=dict(color="orange"),
                            name="break depth (proper)")
    
    trace10 = go.Scatter(x=np.arange(data["contig_length"]),
                            y=-depth_break_all,
                            mode='lines',
                        line=dict(color="orchid"),
                            name="break depth (all)")
    
    ## Start/end count plot
    #count_threshold = 3   # only on the points with > count_threshold start [end] count, plot start [end] data
    array_n_break_start_all = data["array_n_break_start_proper"] + data["array_n_break_start_clipped"]   # TODO: proper>0でall=0のことがある。aggregationの後に再度合計？
    x_break_start_proper, x_break_start_all = find_spike(data["array_n_break_start_proper"], array_n_break_start_all, count_threshold, window_size)
    x_break_start = sorted(list(set(x_break_start_proper) | set(x_break_start_all)))
    array_n_break_end_all = data["array_n_break_end_proper"] + data["array_n_break_end_clipped"]   # TODO: ここら辺色々と怪しいので見直す
    x_break_end_proper, x_break_end_all = find_spike(data["array_n_break_end_proper"], array_n_break_end_all, count_threshold, window_size)
    x_break_end = sorted(list(set(x_break_end_proper) | set(x_break_end_all)))    
    
    trace11 = go.Scatter(x=x_break_start,
                    y=-data["array_n_break_start_proper"][x_break_start],
                        mode="markers",
                        marker=dict(color="darkorange", symbol=8, size=14),
                   name="break start (proper)")
    
    trace12 = go.Scatter(x=x_break_start,
                    y=-array_n_break_start_all[x_break_start],
                        mode="markers",
                        marker=dict(color="purple", symbol=8, size=14),
                   name="break start (all)")
    
    trace13 = go.Scatter(x=x_break_end,
                    y=-data["array_n_break_end_proper"][x_break_end],
                        mode="markers",
                        marker=dict(color="darkorange", symbol=7, size=14),
                   name="break end (proper)")
    
    trace14 = go.Scatter(x=x_break_end,
                    y=-array_n_break_end_all[x_break_end],
                        mode="markers",
                        marker=dict(color="purple", symbol=7, size=14),
                   name="break end (all)")
    
    # Repeat/break annotation
    repeats, breaks = annotate_regions(counts, name)
    y_pos = np.max(depth_all) * 1.1

    shapes = []
    for s, t in repeats:
        shapes.append(dict(type="line", x0=s, y0=y_pos, x1=t, y1=y_pos, line=dict(color="olive")))
    
    for s, t in breaks:
        shapes.append(dict(type="line", x0=s, y0=y_pos, x1=t, y1=y_pos, line=dict(color="red")))
    
    layout = go.Layout(width=width, height=height,
                       margin=go.Margin(l=30, r=30, t=50, b=100, pad=-30),
                           hovermode='closest',
                       title="%s (%d bp window)" % (contig_header, window_size),
                           xaxis=dict(showgrid=False, zeroline=False),
                           yaxis=dict(showgrid=False),
                       shapes=shapes,
                      legend=dict(orientation="h"))
    
    
    figure = go.Figure(data=[trace5, trace6, trace7, trace8, trace1, trace2, trace3, trace4, trace11, trace12, trace13, trace14, trace9, trace10], layout=layout)
    if out_html is not None:
        py.plot(figure, filename=out_html)
    else:
        py.iplot(figure)
