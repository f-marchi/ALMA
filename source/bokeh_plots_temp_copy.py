
from bokeh.layouts import column
from bokeh.models import (ColumnDataSource, CategoricalColorMapper, HoverTool, Label, Span, GroupFilter,
                          CDSView, Legend, LegendItem, Tabs, TabPanel, CustomJS, FactorRange)
from bokeh.plotting import figure, show, output_file
import numpy as np
import pandas as pd

def get_custom_color_palette():
    return [
        '#ff7f0e', '#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#7f7f7f', '#e377c2',
        '#e7ba52', '#bcbd22', '#17becf', '#393b79', '#8c564b', '#f7b6d2', '#c49c94',
        '#a2769e', '#dbdb8d', '#9edae5', '#c5b0d5', '#c7c7c7', '#ff9896', '#637939',
        '#aec7e8', '#ffbb78', '#98df8a', '#7c231e', '#3d6a3d', '#f96502', '#6d3f7d',
        '#6b4423', '#d956a6']

def plot_alma(df, 
            xaxis="PaCMAP 1 of 2", yaxis="PaCMAP 2 of 2",
            x_range=(-45, 40), y_range=(-50, 45), 
            cols=[
                'WHO 2022 Diagnosis', 'AL Epigenomic Subtype', 'Vital Status',
                'AML Epigenomic Risk', 'Risk Group AAML1831', 'Clinical Trial',
                'Race or ethnic group', 'Age (group years)'
                    ],
            save_html=False):
    
    source = ColumnDataSource(df)
    width = 1000
    font_size = "8pt"
    x = 'P(Death) at 5y'
    y = 'Percentile'
    threshold = 0.5

    custom_color_palette = get_custom_color_palette()
    tabs = []

    for col in cols:
        factors = [str(val) for val in df[col].unique() if pd.notnull(val)]
        color_mapper = CategoricalColorMapper(factors=factors, palette=custom_color_palette)

        p2 = figure(title='Acute Leukemia Methylome Atlas', width=width, height=600,
                    tools="pan,wheel_zoom,box_select,reset,save", tooltips=[(str(col), '@{' + str(col) + '}')], 
                    x_axis_label=xaxis, y_axis_label=yaxis,
                    active_drag="box_select", x_range=x_range, y_range=y_range)

        p2.toolbar.logo = None
        p2.toolbar_location = 'above'

        def set_axis_properties(plot, font_size):
            plot.xaxis.axis_label_text_font_size = font_size
            plot.yaxis.axis_label_text_font_size = font_size
            plot.xaxis.axis_label_text_font_style = "normal"
            plot.yaxis.axis_label_text_font_style = "normal"

        set_axis_properties(p2, font_size)

        # Create scatter plot for each factor
        for factor in factors:
            view = CDSView(filter=GroupFilter(column_name=col, group=factor))
            p2.scatter(x=xaxis, y=yaxis, source=source, view=view, 
                       color={'field': col, 'transform': color_mapper}, size=3, alpha=0.8, radius=0.2)

        # Create a list of legend items
        legend_items = [LegendItem(label=factor, renderers=[r]) for factor, r in zip(factors, p2.renderers)]

        # Create a legend
        legend = Legend(items=legend_items, location="top", click_policy="hide",
                        label_text_font_size=font_size, label_text_font_style="normal",
                        glyph_height=15, glyph_width=15, spacing=1)

        # Add the legend to the plot
        p2.add_layout(legend, 'right')

        tab = TabPanel(child=p2, title=col)
        tabs.append(tab)

    alma = Tabs(tabs=tabs, tabs_location='above')

    p1 = figure(title='AML Epigenomic Risk', width=width, height=300,
                tools="xbox_select,reset,save", active_drag='xbox_select',
                x_axis_label=x, y_axis_label="Patient Percentile")
    p1.toolbar.logo = None

    # Add background color to plot1
    p1.quad(left=0.2, right=threshold, bottom=0, top=1, color="#1f77b4", level="underlay", alpha=0.2)
    p1.quad(left=threshold, right=0.8, bottom=0, top=1, color="#ff7f0e", level="underlay", alpha=0.2)

    # Add a label to the line
    label1 = Label(y=0.05, x=threshold + 0.01, text='High Risk', text_font_size='8pt',
                   text_color='#ff7f0e', text_alpha=0.8, text_align='left')
    label2 = Label(y=0.05, x=threshold - 0.01, text='Low Risk', text_font_size='8pt',
                   text_color='#1f77b4', text_alpha=0.8, text_align='right')
    
    p1.add_layout(label1)
    p1.add_layout(label2)

    scatter1 = p1.circle(x, y, source=source, color="steelblue", alpha=0.1, 
                         size=7, hover_alpha=0.5, line_color=None, hover_fill_color="midnightblue",
                         hover_line_color="white", selection_color="midnightblue", selection_alpha=0.7,
                         selection_line_color="white")

    scatter1_hover_tool = HoverTool(renderers=[scatter1], mode='vline', tooltips=None)
    p1.add_tools(scatter1_hover_tool)

    set_axis_properties(p1, font_size)

    # Create the histogram plot for P(Death)
    p3 = figure(title='AML Epigenomic Risk - Select clusters on the map and their prognosis will appear here', width=width, height=300,
                x_axis_label=x, y_axis_label='Frequency', tools="save",)
    p3.toolbar.logo = None

    hist, edges = np.histogram(df[x], bins=100, range=[0.2, 0.8])
    hist_source = ColumnDataSource(data=dict(top=hist, left=edges[:-1], right=edges[1:]))

    hist_renderer = p3.quad(top='top', bottom=0, left='left', right='right', source=hist_source, fill_color="navy", line_color="white", alpha=0.5)

    avg_p_death = df[x].mean()
    avg_line = Span(location=avg_p_death, dimension='height', line_color='black', line_dash='dashed', line_alpha=0.8)
    avg_label = Label(x=avg_p_death, y=max(hist), text=f'Avg: {avg_p_death:.2f}', text_font_size='8pt', text_color='black', text_align='center')
    
    p3.add_layout(avg_line)
    p3.add_layout(avg_label)

    # CustomJS callback to update histogram and average P(Death) based on selection
    callback = CustomJS(args=dict(source=source, hist_source=hist_source, edges=edges, p3=p3, avg_line=avg_line, avg_label=avg_label), code="""
        const indices = source.selected.indices;
        const data = source.data;
        const hist_data = hist_source.data;

        const hist = new Array(edges.length - 1).fill(0);
        let sum_p_death = 0;

        for (let i = 0; i < indices.length; i++) {
            const idx = indices[i];
            const value = data['P(Death) at 5y'][idx];
            sum_p_death += value;
            for (let j = 0; j < edges.length - 1; j++) {
                if (value >= edges[j] && value < edges[j + 1]) {
                    hist[j] += 1;
                    break;
                }
            }
        }

        hist_data['top'] = hist;
        hist_source.change.emit();

        const avg_p_death = sum_p_death / indices.length;
        avg_line.location = avg_p_death;
        avg_label.x = avg_p_death;
        avg_label.text = `Avg: ${avg_p_death.toFixed(2)}`;
        p3.request_render();
    """)

    source.selected.js_on_change('indices', callback)

    if save_html:
        output_file("../data/AML_Epigenomic_Risk.html")


    return show(column(alma, p3, p1))
    

