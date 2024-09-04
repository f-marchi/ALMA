import numpy as np
import pandas as pd
from bokeh.layouts import column
from bokeh.models import (ColumnDataSource, CategoricalColorMapper, Label, Span, GroupFilter,
                          CDSView, Legend, LegendItem, Tabs, TabPanel, CustomJS)
from bokeh.plotting import figure, show, output_file

def get_custom_color_palette():
    return [
        '#ff7f0e', '#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#7f7f7f', '#e377c2',
        '#e7ba52', '#bcbd22', '#17becf', '#393b79', '#8c564b', '#f7b6d2', '#c49c94',
        '#a2769e', '#dbdb8d', '#9edae5', '#c5b0d5', '#c7c7c7', '#ff9896', '#637939',
        '#aec7e8', '#ffbb78', '#98df8a', '#7c231e', '#3d6a3d', '#f96502', '#6d3f7d',
        '#6b4423', '#d956a6']

def set_axis_properties(plot, font_size):
    plot.xaxis.axis_label_text_font_size = font_size
    plot.yaxis.axis_label_text_font_size = font_size
    plot.xaxis.axis_label_text_font_style = "normal"
    plot.yaxis.axis_label_text_font_style = "normal"

def create_scatter_plot(df, col, xaxis, yaxis, x_range, y_range, width, font_size):
    source = ColumnDataSource(df)
    factors = [str(val) for val in df[col].unique() if pd.notnull(val)]
    color_mapper = CategoricalColorMapper(factors=factors, palette=get_custom_color_palette())

    p = figure(title='Acute Leukemia Methylome Atlas', width=width, height=675,
               tools="pan,wheel_zoom,box_select,reset,save", tooltips=[(str(col), f'@{{{col}}}')],
               x_axis_label=xaxis, y_axis_label=yaxis,
               active_drag="box_select", x_range=x_range, y_range=y_range)

    p.toolbar.logo = None
    p.toolbar_location = 'above'
    set_axis_properties(p, font_size)

    for factor in factors:
        view = CDSView(filter=GroupFilter(column_name=col, group=factor))
        p.scatter(x=xaxis, y=yaxis, source=source, view=view,
                  color={'field': col, 'transform': color_mapper}, size=3, alpha=0.8, radius=0.2)

    legend_items = [LegendItem(label=factor, renderers=[r]) for factor, r in zip(factors, p.renderers)]
    legend = Legend(items=legend_items, location="top", click_policy="hide",
                    label_text_font_size=font_size, label_text_font_style="normal",
                    glyph_height=15, glyph_width=15, spacing=1)
    p.add_layout(legend, 'right')

    return p, source
def create_risk_histogram_plot(df, width, font_size, x, y, threshold):
    source = ColumnDataSource(df)
    
    # Create the main figure
    p = figure(title='AML Epigenomic Risk', width=width, height=400,
               tools="box_select,reset,save", active_drag='box_select',
               x_axis_label=x, y_axis_label="Normalized Frequency")
    p.toolbar.logo = None
    set_axis_properties(p, font_size)

    # Add risk zones
    p.quad(left=0.2, right=threshold, bottom=0, top=1, color="#1f77b4", level="underlay", alpha=0.2)
    p.quad(left=threshold, right=0.8, bottom=0, top=1, color="#ff7f0e", level="underlay", alpha=0.2)

    # Add labels
    p.add_layout(Label(y=0.05, x=threshold + 0.01, text='High Risk', text_font_size='8pt',
                       text_color='#ff7f0e', text_alpha=0.8, text_align='left'))
    p.add_layout(Label(y=0.05, x=threshold - 0.01, text='Low Risk', text_font_size='8pt',
                       text_color='#1f77b4', text_alpha=0.8, text_align='right'))

    # Add histogram
    hist, edges = np.histogram(df[x], bins=100, range=[0.2, 0.8])
    hist_normalized = hist / np.max(hist)  # Normalize histogram
    hist_source = ColumnDataSource(data=dict(top=hist_normalized, left=edges[:-1], right=edges[1:]))

    p.quad(top='top', bottom=0, left='left', right='right', source=hist_source,
           fill_color="navy", line_color="white", alpha=0.5)

    # Add average line
    avg_p_death = df[x].mean()
    avg_line = Span(location=avg_p_death, dimension='height', line_color='black', line_dash='dashed', line_alpha=0.8)
    avg_label = Label(x=avg_p_death, y=1, text=f'Avg: {avg_p_death:.2f}',
                      text_font_size='8pt', text_color='black', text_align='center')

    p.add_layout(avg_line)
    p.add_layout(avg_label)

    return p, source, hist_source, avg_line, avg_label, edges

def plot_alma(df, xaxis="PaCMAP 1 of 2", yaxis="PaCMAP 2 of 2",
              x_range=(-45, 40), y_range=(-50, 45),
              cols=['WHO 2022 Diagnosis', 'AL Epigenomic Subtype', 'Vital Status',
                    'AML Epigenomic Risk', 'Risk Group AAML1831', 'Clinical Trial',
                    'Race or ethnic group', 'Age (group years)'],
              save_html=False):

    width = 1000
    font_size = "8pt"
    x = 'P(Death) at 5y'
    y = 'Percentile'
    threshold = 0.5

    tabs = []
    scatter_sources = []
    for col in cols:
        scatter_plot, scatter_source = create_scatter_plot(df, col, xaxis, yaxis, x_range, y_range, width, font_size)
        tabs.append(TabPanel(child=scatter_plot, title=col))
        scatter_sources.append(scatter_source)
    
    alma = Tabs(tabs=tabs, tabs_location='above')

    risk_hist_plot, risk_source, hist_source, avg_line, avg_label, edges = create_risk_histogram_plot(df, width, font_size, x, y, threshold)

    callback = CustomJS(args=dict(scatter_sources=scatter_sources, risk_source=risk_source, 
                                  hist_source=hist_source, edges=edges, avg_line=avg_line, avg_label=avg_label), code="""
        const x_field = 'P(Death) at 5y';
        let all_indices = new Set();

        // Collect indices from all scatter plots and risk plot
        scatter_sources.forEach(source => {
            source.selected.indices.forEach(index => all_indices.add(index));
        });
        risk_source.selected.indices.forEach(index => all_indices.add(index));

        // If no selection, use all indices
        if (all_indices.size === 0) {
            for (let i = 0; i < risk_source.get_length(); i++) {
                all_indices.add(i);
            }
        }

        const hist = new Array(edges.length - 1).fill(0);
        let sum_p_death = 0;
        let count = 0;

        all_indices.forEach(idx => {
            const value = risk_source.data[x_field][idx];
            sum_p_death += value;
            count += 1;
            for (let j = 0; j < edges.length - 1; j++) {
                if (value >= edges[j] && value < edges[j + 1]) {
                    hist[j] += 1;
                    break;
                }
            }
        });

        // Normalize histogram
        const max_hist = Math.max(...hist);
        const hist_normalized = hist.map(h => h / max_hist);

        hist_source.data['top'] = hist_normalized;
        hist_source.change.emit();

        if (count > 0) {
            const avg_p_death = sum_p_death / count;
            avg_line.location = avg_p_death;
            avg_label.x = avg_p_death;
            avg_label.y = 1;  // Set to top of normalized range
            avg_label.text = `Avg: ${avg_p_death.toFixed(2)}`;
        } else {
            avg_line.location = 0;
            avg_label.x = 0;
            avg_label.y = 0;
            avg_label.text = 'No selection';
        }

        // Update the title of the risk_hist_plot
        risk_hist_plot.title.text = count < risk_source.get_length() ? `AML Epigenomic Risk - Selected ${count} patients` : 'AML Epigenomic Risk - All patients';

        // Update scatter plot selections
        scatter_sources.forEach(source => {
            source.selected.indices = Array.from(all_indices);
        });

        // Update risk plot selection
        risk_source.selected.indices = Array.from(all_indices);
    """)

    # Connect the callback to all scatter plots and the risk plot
    for source in scatter_sources:
        source.selected.js_on_change('indices', callback)
    risk_source.selected.js_on_change('indices', callback)

    # Add a reset callback to ensure the histogram is redrawn when reset is clicked
    reset_callback = CustomJS(args=dict(callback=callback), code="""
        callback.execute();
    """)

    for tab in tabs:
        tab.child.js_on_event('reset', reset_callback)
    risk_hist_plot.js_on_event('reset', reset_callback)

    if save_html:
        output_file("../data/AML_Epigenomic_Risk.html")

    return show(column(alma, risk_hist_plot))