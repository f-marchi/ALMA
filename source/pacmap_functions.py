"""
This module implements functions for unsupervised learning algorithms.

"""

import numpy as np
import pandas as pd
import pacmap
from bokeh.plotting import figure, show
from bokeh.transform import factor_cmap
from bokeh.models import Div, Slider, TabPanel, Tabs, Legend, ColumnDataSource
from bokeh.layouts import layout
from bokeh.io import curdoc, output_notebook
from bokeh.embed import json_item

class DataProcessor:
    def __init__(self, train_clinical_data, df_train, clinical_trials, sample_types,
                 cols, n_components, remove_patient_id_duplicates=False,
                 common_prefix=None, df_test=None, test_clinical_data=None):
        
        self.train_clinical_data = train_clinical_data
        self.df_train = df_train
        self.clinical_trials = clinical_trials
        self.sample_types = sample_types
        self.cols = cols
        self.n_components = n_components
        self.remove_patient_id_duplicates = remove_patient_id_duplicates
        self.common_prefix = common_prefix
        self.df_test = df_test
        self.test_clinical_data = test_clinical_data

    def filter_data(self):
        self.df1 = self.train_clinical_data[self.train_clinical_data['Clinical Trial'].isin(self.clinical_trials)]
        self.df2 = self.df1[self.df1['Sample Type'].isin(self.sample_types)]
        if self.remove_patient_id_duplicates:
            self.df3 = self.df2[~self.df2['Patient_ID'].duplicated(keep='last')]
        else:
            self.df3 = self.df2
        self.df_train_filtered = self.df_train[self.df_train.index.isin(self.df3.index)]

    def apply_pacmap(self):

        self.reducer = pacmap.PaCMAP(n_components=self.n_components, n_neighbors=15, MN_ratio=0.4, FP_ratio=16.0, 
                                random_state=42, lr=0.1, num_iters=5000)
        
        embedding = self.reducer.fit_transform(self.df_train_filtered.to_numpy(dtype='float16'))

        if self.common_prefix:
            pacmap.save(self.reducer, self.common_prefix)

        cols = ['PaCMAP ' + str(i+1) for i in range(self.n_components)]

        self.df_embedding = pd.DataFrame(embedding, index=self.df_train_filtered.index, columns=cols)

    def apply_pacmap_test(self):
        if self.df_test is not None:
            embedding_test = self.reducer.transform(self.df_test.to_numpy(dtype='float16'), 
                                                    basis=self.df_train_filtered.to_numpy(dtype='float16').copy())
            cols_test = ['PaCMAP ' + str(i+1) for i in range(self.n_components)]
            self.df_test_embedding = pd.DataFrame(embedding_test, index=self.df_test.index, columns=cols_test)

    def join_labels(self):
        self.df_train = self.df_embedding.join(self.train_clinical_data[self.cols]).reset_index()
        self.df_test = self.df_test_embedding.join(self.test_clinical_data[self.cols]).reset_index()
        self.df = pd.concat([self.df_train, self.df_test]) # Concatenate train and test dataframes


class BokehPlotter:
    from bokeh.themes import Theme

    white_theme = Theme(json={
        "attrs": {
            # "Plot": { "toolbar_location": None },
            # "Grid": { "grid_line_color": None },
            "Axis": {
            #     "axis_line_color": None,
                "major_label_text_color": 'black',
                "major_label_text_font": 'Arial',
            #     "major_tick_line_color": None,
            #     "minor_tick_line_color": None,
            },
            "Legend": {
                "label_text_color": 'black',
                "label_text_font": 'Arial',
            },
            "Title": {
                "text_color": 'black',
                "text_font": 'Arial',
            },
        }
    })

    def __init__(self, df, cols, custom_color_palette, title=None,
                x_range=None, y_range=None, datapoint_size=5, 
                tooltip_dx_cols='WHO 2022 Diagnosis', width=1300, height=800):
        self.df = df
        self.cols = cols
        self.custom_color_palette = custom_color_palette
        self.title = title #+ ', n=' + str(self.df.shape[0])
        self.x_range = x_range or (-50, 50)
        self.y_range = y_range or (-50, 50)
        self.tabs = None
        self.points = None
        self.slider = None
        self.layout = None
        self.datapoint_size = datapoint_size or 5
        self.tooltip_dx_cols = tooltip_dx_cols
        self.width = width
        self.height = height

    def create_figure(self):
        p = figure(title=self.title, 
                width=self.width, height=self.height, sizing_mode="inherit",
                x_axis_label='Longitude (PaCMAP 1)', y_axis_label='Latitude (PaCMAP 2)',
                x_range=self.x_range, y_range=self.y_range,
                tools="pan,wheel_zoom,reset,save", active_drag=None,
                active_scroll="auto",
                tooltips=[("Dx", "@{"+self.tooltip_dx_cols+"}")])
        curdoc().theme = BokehPlotter.white_theme
        return p



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
        self.tabs = Tabs(tabs=[TabPanel(child=self.create_figure(), title=title) for title in self.cols],
                tabs_location='left')

        self.points = [self.create_scatters(tab.child, hue=col) for tab, col in zip(self.tabs.tabs, self.cols)]

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

    def create_json_item(self):
        # This function is used to embed the plot in the web app
        return json_item(self.layout)

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
