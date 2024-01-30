"""
This module implements custom functions for a combined Bokeh plot with two scatter plots and a table linked.

"""  
import numpy as np
import pandas as pd
from bokeh.layouts import column
from bokeh.plotting import figure, show, output_notebook, save

from bokeh.models import (Slider, TabPanel, Tabs, Legend, ColumnDataSource, LegendItem,
                          CDSView, GroupFilter, CategoricalColorMapper, Label, Span,
                          DataTable, HoverTool, TableColumn)

output_notebook()

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

def plot_bokeh(df):

    df2 = df.sort_values(by='P(Dead)').reset_index().reset_index(names=['Percentile'])
    df2['Percentile'] = df2['Percentile'] / len(df2)
    
    source = ColumnDataSource(df2)
    tooltips = [("EpiDx", "@{AL Epigenomic Phenotype}")]
    width = 1000
    tools = "pan,wheel_zoom,box_select,reset,save"
    active_drag = "box_select"
    font_size = "8pt"

    custom_color_palette = get_custom_color_palette()

    cols = ['AL Epigenomic Phenotype', 'Hematopoietic Entity','WHO 2022 Diagnosis', 
            'Vital Status','AML Epigenomic Risk', 'Risk Group AAML1831', 'Clinical Trial' ]

    columns = [TableColumn(field=col, title=col) for col in cols]
    columns.append(TableColumn(field='Gene Fusion', title='Gene Fusion'))
    columns.append(TableColumn(field='Karyotype', title='Karyotype'))


    data_table = DataTable(source=source, columns=columns, editable=True, width=width,
                        index_position=None, height=300)

    p1 = figure(title= 'AML Epigenomic Risk' ,width=width, height=300,
                tools="xbox_select,reset,save", active_drag='xbox_select',
                x_axis_label='AML Epigenomic Risk Score', y_axis_label="est. probability of progression")
    p1.toolbar.logo = None

    threshold = 0.5  # Set your threshold value

    # Add background color to plot1
    p1.quad(left=0, right=threshold, bottom=0, top=1, color="#1f77b4", level="underlay", alpha=0.2)
    p1.quad(left=threshold, right=1, bottom=0, top=1, color="#ff7f0e", level="underlay", alpha=0.2)

    # Add a label to the line
    label1 = Label(y=0.05, x=threshold + 0.01, text='High Risk', text_font_size='8pt',
                   text_color='#ff7f0e', text_alpha=0.8, text_align='left')
    label2 = Label(y=0.05, x=threshold - 0.01, text='Low Risk', text_font_size='8pt',
                   text_color='#1f77b4', text_alpha=0.8, text_align='right')
    p1.add_layout(label1)
    p1.add_layout(label2)


    scatter1 = p1.circle(y='Percentile', x='P(Dead)', source=source, selection_color='#ff7f0e', 
                nonselection_alpha=1.0, color='#1f77b4', size=5, alpha=0.8, hover_color='#ff7f0e',
                hover_alpha=1.0)

    scatter1_hover_tool = HoverTool(renderers=[scatter1], tooltips=tooltips)

    tabs = []

    for col in cols:
        factors = [str(val) for val in df2[col].unique() if pd.notnull(val)]
        color_mapper = CategoricalColorMapper(factors=factors, palette=custom_color_palette)

        p2 = figure(title='Acute Leukemia Methylome Atlas', width=width, height=600,
                    tools="pan,wheel_zoom,box_select,reset,save", tooltips=tooltips, 
                    x_axis_label='Longitude (PaCMAP 1)', y_axis_label='Latitude (PaCMAP 2)',
                    active_drag="box_select", x_range= (-45,40), y_range=(-50,45),)

        p2.toolbar.logo = None
        p2.toolbar_location = 'above'

        def set_axis_properties(plot, font_size):
            plot.xaxis.axis_label_text_font_size = font_size
            plot.yaxis.axis_label_text_font_size = font_size
            plot.xaxis.axis_label_text_font_style = "normal"
            plot.yaxis.axis_label_text_font_style = "normal"

        set_axis_properties(p1, font_size)
        set_axis_properties(p2, font_size)

        # Create scatter plot for each factor
        for factor in factors:
            view = CDSView(filter=GroupFilter(column_name=col, group=factor))
            p2.scatter(x="PaCMAP 1 of 2", y="PaCMAP 2 of 2", source=source, view=view, 
                    color={'field': col, 'transform': color_mapper}, size=3, alpha=0.8, radius=0.2)

        # Create a list of legend items
        legend_items = [LegendItem(label=factor, renderers=[r]) for factor, r in zip(factors, p2.renderers)]

        # Create a legend
        legend = Legend(items=legend_items, location="top", click_policy="hide",
                        label_text_font_size=font_size, label_text_font_style="normal",
                        glyph_height=15, glyph_width=15, spacing=1, margin=5, padding=0,
                        )

        # Add the legend to the plot
        p2.add_layout(legend, 'right')

        tab = TabPanel(child=p2, title=col)
        tabs.append(tab)

    tabs_control = Tabs(tabs=tabs, tabs_location='above')

    p1.add_tools(scatter1_hover_tool)

    return(show(column(tabs_control,p1,data_table)))