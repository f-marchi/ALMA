import numpy as np
import pandas as pd
from bokeh.layouts import column
from bokeh.plotting import figure, show
from bokeh.models import (
    ColumnDataSource, CategoricalColorMapper, HoverTool,
    Label, Span, GroupFilter, CDSView, Legend, LegendItem, Tabs, TabPanel,
    CustomJS, FactorRange
)

def get_custom_color_palette():
    return [
        '#ff7f0e', '#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#7f7f7f',
        '#e377c2', '#e7ba52', '#bcbd22', '#17becf', '#393b79', '#8c564b',
        '#f7b6d2', '#c49c94', '#a2769e', '#dbdb8d', '#9edae5', '#c5b0d5',
        '#c7c7c7', '#ff9896', '#637939', '#aec7e8', '#ffbb78', '#98df8a',
        '#7c231e', '#3d6a3d', '#f96502', '#6d3f7d', '#6b4423', '#d956a6'
    ]

def create_risk_plot(df, source, width, x, threshold, test_sample):
    p = figure(title='AML Epigenomic Risk', width=width, height=300,
               tools="xbox_select,reset,save", active_drag='xbox_select',
               x_axis_label=x, y_axis_label="Patient Percentile")
    p.toolbar.logo = None

    p.quad(left=0.15, right=threshold, bottom=0, top=1, color="#1f77b4", level="underlay", alpha=0.2)
    p.quad(left=threshold, right=0.85, bottom=0, top=1, color="#ff7f0e", level="underlay", alpha=0.2)

    for label_text, x_pos, color, align in [('High Risk', threshold + 0.01, '#ff7f0e', 'left'),
                                            ('Low Risk', threshold - 0.01, '#1f77b4', 'right')]:
        p.add_layout(Label(y=0.05, x=x_pos, text=label_text, text_font_size='8pt',
                           text_color=color, text_alpha=0.8, text_align=align))

    scatter = p.circle(x, 'Percentile', source=source, color="steelblue", alpha=0.1, 
                       size=7, hover_alpha=0.5, line_color=None, hover_fill_color="midnightblue",
                       hover_line_color="white", selection_color="midnightblue", selection_alpha=0.7,
                       selection_line_color="white")

    p.add_tools(HoverTool(renderers=[scatter], mode='vline', tooltips=None))

    if test_sample:
        vline = Span(location=df.loc[test_sample][x], dimension='height', line_color='black', line_dash='dashed', line_alpha=0.8)
        p.renderers.extend([vline])
        p.star(x=df.loc[test_sample][x], y=0, size=15, color="black", alpha=0.9, 
               legend_label=f'{test_sample}, {df.loc[test_sample]["AML Epigenomic Risk"]} Epigenomic Risk ({df.loc[test_sample][x]:.2f})',
               line_color="black", line_width=1)
        p.legend.location = "bottom_right"
        p.legend.click_policy = "hide"

    return p

def create_histogram_plot(df, source, width, x):
    p = figure(title='AML Epigenomic Risk - Select clusters on the map and their prognosis will appear here', 
               width=width, height=300, x_axis_label=x, y_axis_label='Frequency', tools="save,reset")
    p.toolbar.logo = None

    hist, edges = np.histogram(df[x], bins=50, range=[0.15, 0.85])
    hist_source = ColumnDataSource(data=dict(top=hist, left=edges[:-1], right=edges[1:]))
    p.quad(top='top', bottom=0, left='left', right='right', source=hist_source, fill_color="navy", line_color="white", alpha=0.5)

    callback = CustomJS(args=dict(source=source, hist_source=hist_source, edges=edges, p=p), code="""
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
        p.request_render();
    """)

    source.selected.js_on_change('indices', callback)

    return p

def create_race_risk_plot(df, source, width):
    race_risk_counts = df.groupby(['Race or ethnic group', 'AML Epigenomic Risk']).size().unstack().fillna(0)
    race_totals = race_risk_counts.sum(axis=1)
    race_risk_counts = race_risk_counts.div(race_totals, axis=0) * 100
    race_categories = list(race_risk_counts.index)
    race_hist_source = ColumnDataSource(data=dict(race=race_categories, high=race_risk_counts['High'], low=race_risk_counts['Low'], count=race_totals))

    p = figure(title='Race or Ethnic Group by AML Epigenomic Risk', width=width, height=300,
               x_range=FactorRange(*race_categories), y_axis_label='Percentage', y_range=(0, 100), tools="save, reset")
    p.toolbar.logo = None

    hover = HoverTool(tooltips=[("High Risk", "@high{0.0}%"), ("Low Risk", "@low{0.0}%"), ("Count", "@count")])
    p.add_tools(hover)

    p.vbar_stack(['low', 'high'], x='race', width=0.7, color=["#1f77b4", "#ff7f0e"], source=race_hist_source,
                 legend_label=['Low AML Epigenomic Risk', 'High AML Epigenomic Risk'], line_color="white", alpha=0.5)

    p.legend.location = "top_left"

    callback = CustomJS(args=dict(source=source, race_hist_source=race_hist_source, race_categories=race_categories), code="""
        const indices = source.selected.indices;
        const data = source.data;
        const race_hist_data = race_hist_source.data;

        const race_risk_counts = {};
        for (let i = 0; i < race_categories.length; i++) {
            race_risk_counts[race_categories[i]] = {High: 0, Low: 0};
        }

        for (let i = 0; i < indices.length; i++) {
            const idx = indices[i];
            const race = data['Race or ethnic group'][idx];
            const risk = data['AML Epigenomic Risk'][idx];
            if (race in race_risk_counts) {
                race_risk_counts[race][risk] += 1;
            }
        }

        for (let i = 0; i < race_categories.length; i++) {
            const race = race_categories[i];
            const high = race_risk_counts[race]['High'];
            const low = race_risk_counts[race]['Low'];
            const total = high + low;
            race_hist_data['high'][i] = (total > 0) ? (high / total) * 100 : 0;
            race_hist_data['low'][i] = (total > 0) ? (low / total) * 100 : 0;
            race_hist_data['count'][i] = total;
        }

        race_hist_source.change.emit();
    """)

    source.selected.js_on_change('indices', callback)

    return p

def create_scatter_plot(df, source, col, x_range, y_range, xaxis, yaxis, test_sample):
    factors = [str(val) for val in df[col].unique() if pd.notnull(val)]
    color_mapper = CategoricalColorMapper(factors=factors, palette=get_custom_color_palette())

    p = figure(title='Acute Leukemia Methylome Atlas', width=1000, height=675,
               tools="pan,wheel_zoom,box_select,reset,save", tooltips=[(str(col), '@{' + str(col) + '}')], 
               x_axis_label=xaxis, y_axis_label=yaxis,
               active_drag="box_select", x_range=x_range, y_range=y_range)

    p.toolbar.logo = None
    p.toolbar_location = 'above'

    for axis in [p.xaxis, p.yaxis]:
        axis.axis_label_text_font_size = "8pt"
        axis.axis_label_text_font_style = "normal"

    for factor in factors:
        view = CDSView(filter=GroupFilter(column_name=col, group=factor))
        p.scatter(x=xaxis, y=yaxis, source=source, view=view, 
                  color={'field': col, 'transform': color_mapper}, size=3, alpha=0.8, radius=0.2)

    if test_sample:
        for dim, loc in [('height', df.loc[test_sample][xaxis]), ('width', df.loc[test_sample][yaxis])]:
            p.renderers.extend([Span(location=loc, dimension=dim, line_color="black", line_dash='dashed', line_alpha=0.8)])
        p.star(x=df.loc[test_sample][xaxis], y=df.loc[test_sample][yaxis],
               size=15, color="black", alpha=0.9, legend_label=f'Sample: {test_sample}\nPrediction: {df.loc[test_sample]["AL Epigenomic Subtype"]}',
               line_color="black", line_width=1)
        p.legend.click_policy = "hide"

    legend = Legend(items=[LegendItem(label=factor, renderers=[r]) for factor, r in zip(factors, p.renderers)],
                    location="top", click_policy="hide",
                    label_text_font_size="8pt", label_text_font_style="normal",
                    glyph_height=15, glyph_width=15, spacing=1)

    p.add_layout(legend, 'right')

    return p

def plot_alma(df, test_sample=None, 
              xaxis="PaCMAP 1 of 2", yaxis="PaCMAP 2 of 2",
              x_range=(-45, 45), y_range=(-45, 45), 
              cols=['AL Epigenomic Subtype','WHO 2022 Diagnosis','Hematopoietic Entity', 
                    'Vital Status', 'AML Epigenomic Risk', 
                    'Risk Group AAML1831', 'Clinical Trial',
                    'Race or ethnic group', 'Age (group years)'],
              save_html=False):
    
    source = ColumnDataSource(df)
    width = 1000
    x = 'P(Death) at 5y'
    threshold = 0.5

    risk_plot = create_risk_plot(df, source, width, x, threshold, test_sample)
    histogram_plot = create_histogram_plot(df, source, width, x)
    race_risk_plot = create_race_risk_plot(df, source, width)

    tabs = []
    for col in cols:
        scatter_plot = create_scatter_plot(df, source, col, x_range, y_range, xaxis, yaxis, test_sample)
        tabs.append(TabPanel(child=scatter_plot, title=col))

    tabs_control = Tabs(tabs=tabs, tabs_location='above')

    if save_html:
        from bokeh.plotting import output_file
        output_file("../data/ALMA.html")

    layout = column(tabs_control, histogram_plot, risk_plot, race_risk_plot)

    return show(layout)