import pickle
import numpy as np
import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go

plt.style.use('ggplot')


class Plotter:
    def __init__(self):
        self.contigs = pickle.load(open("contigs.pkl", 'rb'))
        self.reads = pickle.load(open("reads.pkl", 'rb'))
        self.mappings = pickle.load(open("mappings.pkl", 'rb'))
        self.counts = pickle.load(open("counts.pkl", 'rb'))
        self.depths = pickle.load(open("depths.pkl", 'rb'))
        self.spikes = pickle.load(open("spikes.pkl", 'rb'))
        self.annotations = pickle.load(open("annotations.pkl", 'rb'))

    def plot(self, contig_id, width=1000, height=600, out_html=None):
        contig_length = self.contigs.loc[contig_id, "length"]

        x_start_proper = self.spikes.loc[(contig_id, "proper"), "array_n_start"]
        x_start_clipped = self.spikes.loc[(contig_id, "clipped"), "array_n_start"]
        x_end_proper = self.spikes.loc[(contig_id, "proper"), "array_n_end"]
        x_end_clipped = self.spikes.loc[(contig_id, "clipped"), "array_n_end"]
        #x_start = sorted(list(set(x_start_proper) | set(x_start_clipped)))
        #x_end = sorted(list(set(x_end_proper) | set(x_end_clipped)))
        x_break_start_proper = self.spikes.loc[(contig_id, "proper"), "array_n_break_start"]
        x_break_start_clipped = self.spikes.loc[(contig_id, "clipped"), "array_n_break_start"]
        x_break_end_proper = self.spikes.loc[(contig_id, "proper"), "array_n_break_end"]
        x_break_end_clipped = self.spikes.loc[(contig_id, "clipped"), "array_n_break_end"]
        #x_break_start = sorted(list(set(x_break_start_proper) | set(x_break_start_clipped)))
        #x_break_end = sorted(list(set(x_break_end_proper) | set(x_break_end_clipped)))
    
        trace1 = go.Scatter(x=[0,
                               contig_length - 1],
                            y=[self.counts.loc[(contig_id, "proper"), "n_contains"],
                               self.counts.loc[(contig_id, "proper"), "n_contains"]],
                            mode="lines",
                            line=dict(color="grey"),
                            name="contains")
        
        trace2 = go.Scatter(x=np.arange(contig_length),
                            y=self.depths.loc[(contig_id, "proper"), "depth"],
                            mode='lines',
                            line=dict(color="dodgerblue"),
                            name="depth (proper)")
    
        trace3 = go.Scatter(x=np.arange(contig_length),
                            y=self.depths.loc[(contig_id, "clipped"), "depth"],
                            mode='lines',
                            line=dict(color="limegreen"),
                            name="depth (clipped)")
        
        trace4 = go.Scatter(x=[0,
                               0,
                               contig_length - 1,
                               contig_length - 1],
                            y=[self.counts.loc[(contig_id, "proper"), "n_dovetail_start"],
                               self.counts.loc[(contig_id, "clipped"), "n_dovetail_start"],
                               self.counts.loc[(contig_id, "proper"), "n_dovetail_end"],
                               self.counts.loc[(contig_id, "clipped"), "n_dovetail_end"]],
                            mode="markers",
                            marker=dict(color=["blue", "green", "blue", "green"], symbol=24, size=14),
                            name="n_dovetail")
        
        trace5 = go.Scatter(x=x_start_proper,
                            y=(self.counts.loc[(contig_id, "proper"), "array_n_start"][x_start_proper]
                               if len(x_start_proper) > 0 else []),
                            mode="markers",
                            marker=dict(color="blue", symbol=8, size=14),
                            name="start (proper)")
        
        trace6 = go.Scatter(x=x_start_clipped,
                            y=(self.counts.loc[(contig_id, "clipped"), "array_n_start"][x_start_clipped]
                               if len(x_start_clipped) > 0 else []),
                            mode="markers",
                            marker=dict(color="green", symbol=8, size=14),
                            name="start (clipped)")
        
        trace7 = go.Scatter(x=x_end_proper,
                            y=(self.counts.loc[(contig_id, "proper"), "array_n_end"][x_end_proper]
                               if len(x_end_proper) > 0 else []),
                            mode="markers",
                            marker=dict(color="blue", symbol=7, size=14),
                            name="end (proper)")
        
        trace8 = go.Scatter(x=x_end_clipped,
                            y=(self.counts.loc[(contig_id, "clipped"), "array_n_end"][x_end_clipped]
                               if len(x_end_clipped) > 0 else []),
                            mode="markers",
                            marker=dict(color="green", symbol=7, size=14),
                            name="end (clipped)")
        
        trace9 = go.Scatter(x=np.arange(contig_length),
                            y=-self.depths.loc[(contig_id, "proper"), "break_depth"],
                            mode='lines',
                            line=dict(color="orange"),
                            name="break depth (proper)")
        
        trace10 = go.Scatter(x=np.arange(contig_length),
                             y=-self.depths.loc[(contig_id, "clipped"), "break_depth"],
                             mode='lines',
                             line=dict(color="orchid"),
                             name="break depth (clipped)")
        
        trace11 = go.Scatter(x=x_break_start_proper,
                             y=(-self.counts.loc[(contig_id, "proper"), "array_n_break_start"][x_break_start_proper]
                                if len(x_break_start_proper) > 0 else []),
                             mode="markers",
                             marker=dict(color="darkorange", symbol=8, size=14),
                             name="break start (proper)")
        
        trace12 = go.Scatter(x=x_break_start_clipped,
                             y=(-self.counts.loc[(contig_id, "clipped"), "array_n_break_start"][x_break_start_clipped]
                                if len(x_break_start_clipped) > 0 else []),
                             mode="markers",
                             marker=dict(color="purple", symbol=8, size=14),
                             name="break start (clipped)")
        
        trace13 = go.Scatter(x=x_break_end_proper,
                             y=(-self.counts.loc[(contig_id, "proper"), "array_n_break_end"][x_break_end_proper]
                                if len(x_break_end_proper) > 0 else []),
                             mode="markers",
                             marker=dict(color="darkorange", symbol=7, size=14),
                             name="break end (proper)")
        
        trace14 = go.Scatter(x=x_break_end_clipped,
                             y=(-self.counts.loc[(contig_id, "clipped"), "array_n_break_end"][x_break_end_clipped]
                                if len(x_break_end_clipped) > 0 else []),
                             mode="markers",
                             marker=dict(color="purple", symbol=7, size=14),
                             name="break end (clipped)")
        
        # Repeat/break annotation
        repeats_proper = self.annotations.loc[(contig_id, "proper"), "repeats"]
        repeats_clipped = self.annotations.loc[(contig_id, "clipped"), "repeats"]
        breaks_proper = self.annotations.loc[(contig_id, "proper"), "breaks"]
        breaks_clipped = self.annotations.loc[(contig_id, "clipped"), "breaks"]
        y_pos = np.max([np.max(self.depths.loc[(contig_id, "proper"), "depth"]),
                        np.max(self.depths.loc[(contig_id, "clipped"), "depth"])]) * 1.1
        
        shapes = []
        for s, t in repeats_proper:
            shapes.append(dict(type="line", x0=s, y0=y_pos, x1=t, y1=y_pos, line=dict(color="black")))
        for s, t in repeats_clipped:
            shapes.append(dict(type="line", x0=s, y0=y_pos, x1=t, y1=y_pos, line=dict(color="olive")))    
        for s, t in breaks_proper:
            shapes.append(dict(type="line", x0=s, y0=y_pos, x1=t, y1=y_pos, line=dict(color="pink")))
        for s, t in breaks_clipped:
            shapes.append(dict(type="line", x0=s, y0=y_pos, x1=t, y1=y_pos, line=dict(color="red")))
            
        layout = go.Layout(width=width, height=height,
                           margin=go.layout.Margin(l=30, r=30, t=50, b=100, pad=0),
                           hovermode='closest',
                           title=f"{self.contigs.loc[contig_id, 'header']}",
                           xaxis=dict(showgrid=False, zeroline=False),
                           yaxis=dict(showgrid=False),
                           shapes=shapes,
                           legend=dict(orientation="h"))
            
        figure = go.Figure(data=[trace5,
                                 trace6,
                                 trace7,
                                 trace8,
                                 trace1,
                                 trace2,
                                 trace3,
                                 trace4,
                                 trace11,
                                 trace12,
                                 trace13,
                                 trace14,
                                 trace9,
                                 trace10],
                           layout=layout)
        
        if out_html is not None:
            py.plot(figure, filename=out_html)
        else:
            py.iplot(figure)
