import numpy as np
import pandas as pd
from bokeh.layouts import column, row
from bokeh.plotting import figure, show
from bokeh.models import (
    ColumnDataSource, CategoricalColorMapper, HoverTool,
    Label, Span, GroupFilter, CDSView, Legend, LegendItem, Tabs, TabPanel,
    CustomJS, FactorRange, Select
)
from bokeh.io import output_file

def get_custom_color_palette():
    return [
        '#ff7f0e', '#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#7f7f7f',
        '#e377c2', '#e7ba52', '#bcbd22', '#17becf', '#393b79', '#8c564b',
        '#f7b6d2', '#c49c94', '#a2769e', '#dbdb8d', '#9edae5', '#c5b0d5',
        '#c7c7c7', '#ff9896', '#637939', '#aec7e8', '#ffbb78', '#98df8a',
        '#7c231e', '#3d6a3d', '#f96502', '#6d3f7d', '#6b4423', '#d956a6'
    ]

def create_risk_plot(source, width, x, y, threshold, test_sample):
    p = figure(title='Risk Plot', width=width, height=300,
               tools="xbox_select,reset,save", active_drag='xbox_select',
               x_axis_label=x, y_axis_label="Patient Percentile")
    p.toolbar.logo = None

    p.quad(left=0.15, right=threshold, bottom=0, top=1, color="#1f77b4", level="underlay", alpha=0.2)
    p.quad(left=threshold, right=0.85, bottom=0, top=1, color="#ff7f0e", level="underlay", alpha=0.2)

    for label_text, x_pos, color, align in [('High Risk', threshold + 0.01, '#ff7f0e', 'left'),
                                            ('Low Risk', threshold - 0.01, '#1f77b4', 'right')]:
        p.add_layout(Label(y=0.05, x=x_pos, text=label_text, text_font_size='8pt',
                           text_color=color, text_alpha=0.8, text_align=align))

    scatter = p.circle(x, y, source=source, color="steelblue", alpha=0.1, 
                       size=7, hover_alpha=0.5, line_color=None, hover_fill_color="midnightblue",
                       hover_line_color="white", selection_color="midnightblue", selection_alpha=0.7,
                       selection_line_color="white")

    p.add_tools(HoverTool(renderers=[scatter], mode='vline', tooltips=None))

    if test_sample:
        vline = Span(location=source.data[x][test_sample], dimension='height', line_color='black', line_dash='dashed', line_alpha=0.8)
        p.renderers.extend([vline])
        p.star(x=source.data[x][test_sample], y=0, size=15, color="black", alpha=0.9, 
               legend_label=f'{test_sample}, {source.data[y][test_sample]} Risk ({source.data[x][test_sample]:.2f})',
               line_color="black", line_width=1)
        p.legend.location = "bottom_right"
        p.legend.click_policy = "hide"

    return p

def create_histogram_plot(source, width, x):
    p = figure(title='Risk Distribution - Select clusters on the map and their risk will appear here', 
               width=width, height=300, x_axis_label=x, y_axis_label='Frequency', tools="save")
    p.toolbar.logo = None

    hist, edges = np.histogram(source.data[x], bins=50, range=[0.15, 0.85])
    hist_source = ColumnDataSource(data=dict(top=hist, left=edges[:-1], right=edges[1:]))
    p.quad(top='top', bottom=0, left='left', right='right', source=hist_source, fill_color="navy", line_color="white", alpha=0.5)

    return p, hist_source

def create_category_risk_plot(source, width, category, risk):
    category_risk_counts = pd.DataFrame(source.data).groupby([category, risk]).size().unstack().fillna(0)
    category_totals = category_risk_counts.sum(axis=1)
    category_risk_counts = category_risk_counts.div(category_totals, axis=0) * 100
    categories = list(category_risk_counts.index)
    category_hist_source = ColumnDataSource(data=dict(category=categories, high=category_risk_counts['High'], low=category_risk_counts['Low'], count=category_totals))

    p = figure(title=f'{category} by Risk', width=width, height=300,
               x_range=FactorRange(*categories), y_axis_label='Percentage', y_range=(0, 100), tools="save")
    p.toolbar.logo = None

    hover = HoverTool(tooltips=[("High Risk", "@high{0.0}%"), ("Low Risk", "@low{0.0}%"), ("Count", "@count")])
    p.add_tools(hover)

    p.vbar_stack(['low', 'high'], x='category', width=0.7, color=["#1f77b4", "#ff7f0e"], source=category_hist_source,
                 legend_label=['Low Risk', 'High Risk'], line_color="white", alpha=0.5)

    p.legend.location = "top_left"

    return p, category_hist_source

def create_scatter_plot(source, col, x_range, y_range, xaxis, yaxis, test_sample):
    factors = [str(val) for val in set(source.data[col]) if pd.notnull(val)]
    color_mapper = CategoricalColorMapper(factors=factors, palette=get_custom_color_palette())

    p = figure(title='Acute Leukemia Methylome Atlas', width=1000, height=600,
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
        p.scatter(xaxis, yaxis, source=source, view=view, 
                  color={'field': col, 'transform': color_mapper}, size=3, alpha=0.8, radius=0.2)

    if test_sample is not None:
        for dim, loc in [('height', source.data[xaxis][test_sample]), ('width', source.data[yaxis][test_sample])]:
            p.renderers.extend([Span(location=loc, dimension=dim, line_color="black", line_dash='dashed', line_alpha=0.8)])
        p.star(x=source.data[xaxis][test_sample], y=source.data[yaxis][test_sample],
               size=15, color="black", alpha=0.9, legend_label=f'Sample: {test_sample}\nPrediction: {source.data[col][test_sample]}',
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
                    'Vital Status', 'AML Epigenomic Risk', 'MethylScoreAML Categorical',
                    'Risk Group AAML1831', 'Clinical Trial',
                    'Race or ethnic group', 'Age (group years)'],
              save_html=True):
    
    source = ColumnDataSource(df)
    width = 1000
    threshold = 0.5

    risk_options = [('AML Epigenomic Risk', 'P(Death) at 5y'), ('MethylScoreAML Categorical', 'MethylScoreAML')]
    risk_select = Select(title="Select Risk Metric:", value=risk_options[0][0], options=[opt[0] for opt in risk_options])

    risk_plot = create_risk_plot(source, width, risk_options[0][1], 'Percentile', threshold, test_sample)
    histogram_plot, hist_source = create_histogram_plot(source, width, risk_options[0][1])
    category_risk_plot, category_hist_source = create_category_risk_plot(source, width, 'Race or ethnic group', risk_options[0][0])

    callback = CustomJS(args=dict(source=source, risk_plot=risk_plot, histogram_plot=histogram_plot, 
                                  category_risk_plot=category_risk_plot, risk_options=risk_options,
                                  hist_source=hist_source, category_hist_source=category_hist_source,
                                  threshold=threshold), code="""
        const selected_risk = cb_obj.value;
        let x, y;
        for (let option of risk_options) {
            if (option[0] === selected_risk) {
                x = option[1];
                y = option[0];
                break;
            }
        }
        
        // Update risk plot
        risk_plot.xaxis[0].axis_label = x;
        risk_plot.title.text = 'Risk Plot: ' + y;
        
        // Update histogram plot
        histogram_plot.xaxis[0].axis_label = x;
        histogram_plot.title.text = y + ' Distribution - Select clusters on the map and their risk will appear here';
        
        // Update category risk plot
        category_risk_plot.title.text = 'Race or ethnic group by ' + y;
        
        // Update histogram data
        const hist = new Array(hist_source.data.left.length).fill(0);
        const data = source.data[x];
        for (let i = 0; i < data.length; i++) {
            for (let j = 0; j < hist_source.data.left.length; j++) {
                if (data[i] >= hist_source.data.left[j] && data[i] < hist_source.data.right[j]) {
                    hist[j]++;
                    break;
                }
            }
        }
        hist_source.data['top'] = hist;
        
        // Update category risk data
        const categories = category_hist_source.data.category;
        const high = new Array(categories.length).fill(0);
        const low = new Array(categories.length).fill(0);
        const count = new Array(categories.length).fill(0);
        for (let i = 0; i < source.data[y].length; i++) {
            const cat_index = categories.indexOf(source.data['Race or ethnic group'][i]);
            if (cat_index !== -1) {
                if (source.data[y][i] === 'High') {
                    high[cat_index]++;
                } else {
                    low[cat_index]++;
                }
                count[cat_index]++;
            }
        }
        for (let i = 0; i < categories.length; i++) {
            category_hist_source.data.high[i] = count[i] > 0 ? (high[i] / count[i]) * 100 : 0;
            category_hist_source.data.low[i] = count[i] > 0 ? (low[i] / count[i]) * 100 : 0;
            category_hist_source.data.count[i] = count[i];
        }
        
        // Trigger updates
        hist_source.change.emit();
        category_hist_source.change.emit();
    """)

    risk_select.js_on_change('value', callback)

    tabs = []
    for col in cols:
        scatter_plot = create_scatter_plot(source, col, x_range, y_range, xaxis, yaxis, test_sample)
        tabs.append(TabPanel(child=scatter_plot, title=col))

    tabs_control = Tabs(tabs=tabs, tabs_location='above')

    layout = column(risk_select, tabs_control, histogram_plot, risk_plot, category_risk_plot)

    if save_html:
        output_file("ALMA.html")

    return show(layout)