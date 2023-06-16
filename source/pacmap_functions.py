"""
This module implements functions for unsupervised learning algorithms.

"""

import numpy as np
import pandas as pd
import pacmap
import warnings
from bokeh.plotting import figure, show
from bokeh.transform import factor_cmap
from bokeh.models import Div, Slider, TabPanel, Tabs, Legend, ColumnDataSource
from bokeh.layouts import layout
from bokeh.io import curdoc, output_notebook

class DataProcessor:
    def __init__(self, df_labels, df_methyl, clinical_trials, sample_types, cols, remove_duplicates=False):
        self.df_labels = df_labels
        self.df_methyl = df_methyl
        self.clinical_trials = clinical_trials
        self.sample_types = sample_types
        self.cols = cols
        self.remove_duplicates = remove_duplicates

    def filter_data(self):
        self.df1 = self.df_labels[self.df_labels['Clinical Trial'].isin(self.clinical_trials)]
        self.df2 = self.df1[self.df1['Sample Type'].isin(self.sample_types)]
        if self.remove_duplicates:
            self.df3 = self.df2[~self.df2['Patient_ID'].duplicated(keep='last')]
        else:
            self.df3 = self.df2
        self.df_methyl_filtered = self.df_methyl[self.df_methyl.index.isin(self.df3.index)].iloc[:, 1:]

    def apply_pacmap(self):

        # Define a filter to ignore the specific UserWarning
        warnings.filterwarnings("ignore", category=UserWarning, message="Warning: random state is set to")

        reducer = pacmap.PaCMAP(n_components=2, n_neighbors=15, MN_ratio=0.4, FP_ratio=16.0, 
                                random_state=42, lr=0.1, num_iters=5000)
        embedding = reducer.fit_transform(self.df_methyl_filtered.to_numpy(dtype='float16'))

        self.df_embedding = pd.DataFrame(embedding, index=self.df_methyl_filtered.index, columns=['PaCMAP 1', 'PaCMAP 2'])

    def join_labels(self):
        self.df = self.df_embedding.join(self.df_labels[self.cols]).reset_index()


class BokehPlotter:
    def __init__(self, df, cols, custom_color_palette, title=None, x_range=None, y_range=None, datapoint_size=5):
        self.df = df
        self.cols = cols
        self.custom_color_palette = custom_color_palette
        self.title = title + ', n=' + str(self.df.shape[0])
        self.x_range = x_range or (-50, 50)
        self.y_range = y_range or (-50, 50)
        self.tabs = None
        self.points = None
        self.slider = None
        self.layout = None
        self.datapoint_size = datapoint_size or 5

    def create_figure(self):
        return figure(title=self.title,
                      width=1200, height=600, sizing_mode='fixed',
                      x_axis_label='Longitude (PaCMAP 1)', y_axis_label='Latitude (PaCMAP 2)',
                      x_range=self.x_range, y_range=self.y_range,
                      tools="pan,wheel_zoom,reset,save", active_drag="pan",
                      active_scroll="auto",
                      tooltips=[("Karyotype", "@Karyotype"),
                                ("Gene Fusion", "@{Gene Fusion}"),
                                ("Patient ID", "@{Patient_ID}")])

    def create_scatters(self, p, hue):
        df = self.df[~self.df[hue].isna()]  # Filter out rows with NaN values for the hue column
        filtered_dfs = [df[df[hue] == val] for val in df[hue].value_counts().sort_values(ascending=False).index.to_list()]
        
        renderers = []
        items = []
        for i in range(len(filtered_dfs)):
            name = filtered_dfs[i][hue].head(1).values[0]
            color = self.custom_color_palette[i % len(self.custom_color_palette)]
            source = ColumnDataSource(filtered_dfs[i])
            r = p.scatter(x="PaCMAP 1", y="PaCMAP 2", source=source,
                         fill_alpha=0.8, size=self.datapoint_size,
                         color=color)
            renderers.append(r)
            items.append((name, [r]))

        return renderers, items

    def create_tabs_and_points(self):
        self.tabs = Tabs(tabs=[TabPanel(child=self.create_figure(), title=title) for title in self.cols[:-3]],
                tabs_location='left')

        self.points = [self.create_scatters(tab.child, hue=col) for tab, col in zip(self.tabs.tabs, self.cols[:-3])]

    def finalize_tabs(self):
        for p, (renderers, items) in zip(self.tabs.tabs, self.points):
            p.child.toolbar.logo = None
            p.child.toolbar_location = 'above'
            legend = Legend(items=items, location='top_left')
            p.child.add_layout(legend, 'right')
            p.child.legend.click_policy = 'hide'
        for i in range(len(self.tabs.tabs)):
            self.tabs.tabs[i].child.legend.title = self.tabs.tabs[i].title
            self.tabs.tabs[i].child.output_backend = "svg"

    def create_slider(self):
        self.slider = Slider(title="Adjust datapoint size", start=0, end=10, step=1, value=self.points[0][0][0].glyph.size)
        for i in range(len(self.points)): 
            for r in self.points[i][0]: 
                self.slider.js_link("value", r.glyph, "size")

    def create_layout(self):
        div = Div(text="""<br>""", width=1000, height=10)
        self.layout = layout([[[div, self.tabs, self.slider]]])

    def plot(self):
        self.create_tabs_and_points()
        self.finalize_tabs()
        self.create_slider()
        self.create_layout()
        show(self.layout)

def get_custom_color_palette():
    list = [
    '#1f77b4',  # Vivid blue
    '#ff7f0e',  # Vivid orange 
    '#2ca02c',  # Vivid green
    '#d62728',  # Vivid red
    '#9467bd',  # Vivid purple 
    '#7f7f7f',  # Medium gray
    '#e377c2',  # Pink
    '#e7ba52',  # Light orange
    '#bcbd22',  # Olive
    '#17becf',  # Light blue
    '#393b79',  # Dark blue
    '#8c564b',  # Brown
    '#f7b6d2',  # Light pink
    '#c49c94',  # Light brown
    '#a2769e',   # Soft purple
    '#dbdb8d',  # Pale yellow
    '#9edae5',  # Pale cyan
    '#c5b0d5',  # Pale purple
    '#c7c7c7',  # Light gray
    '#ff9896',  # Light red
    '#637939',  # Dark olive
    '#aec7e8',  # Light blue
    '#ffbb78',  # Light orange
    '#98df8a',  # Light green
    '#7c231e',  # Dark red
    '#3d6a3d',  # Dark green
    '#f96502',  # Deep orange
    '#6d3f7d',  # Deep purple
    '#6b4423',  # Dark brown
    '#d956a6'   # Hot pink
    ]
    return list
