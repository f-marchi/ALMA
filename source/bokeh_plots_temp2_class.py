from bokeh.layouts import column
from bokeh.models import (ColumnDataSource, TableColumn, DataTable, HoverTool, Label, Span, 
                          CustomJS, FactorRange, GroupFilter, CDSView, Legend, LegendItem, TabPanel, Tabs, 
                          CategoricalColorMapper)
from bokeh.plotting import figure, show, save, output_notebook
import numpy as np
import pandas as pd

class Linked_ALMA_Plotter:
    def __init__(self, df, xaxis="PaCMAP 1 of 2", yaxis="PaCMAP 2 of 2",
                 x_range=(-45, 40), y_range=(-50, 45), cols=None, custom_color_palette=None):
        if cols is None:
            cols = ['AL Epigenomic Subtype', 'Hematopoietic Entity', 'WHO 2022 Diagnosis', 
                    'Vital Status', 'AML Epigenomic Risk', 'Risk Group AAML1831', 'Clinical Trial',
                    'Race or ethnic group', 'Age (group years)']
        self.df = df
        self.xaxis = xaxis
        self.yaxis = yaxis
        self.x_range = x_range
        self.y_range = y_range
        self.cols = cols
        self.custom_color_palette = custom_color_palette if custom_color_palette else self.get_custom_color_palette()
        self.df2 = self.prepare_data()
        self.source = ColumnDataSource(self.df2)
        self.width = 1000

    def get_custom_color_palette(self):
        return [
    '#ff7f0e',  # Vivid orange
    '#1f77b4',  # Vivid blue 
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
    
    def prepare_data(self):
        df_px = self.df[~self.df['P(Death)'].isna()]
        df_px = df_px.sort_values(by='P(Death)').reset_index().reset_index(names=['Percentile']).set_index('index')
        df_px['Percentile'] = df_px['Percentile'] / len(df_px['Percentile'])
        return self.df.join(df_px[['Percentile']])

    def create_data_table(self):
        columns = [TableColumn(field=col, title=col) for col in self.cols] + [
            TableColumn(field='Gene Fusion', title='Gene Fusion'),
            TableColumn(field='Karyotype', title='Karyotype')
        ]
        return DataTable(source=self.source, columns=columns, editable=True, width=self.width, index_position=None, height=300)

    def create_risk_plot(self, threshold=0.5):
        p = figure(title='AML Epigenomic Risk', width=self.width, height=300,
                   tools="xbox_select,reset,save", active_drag='xbox_select',
                   x_axis_label='P(Death)', y_axis_label="Patient Percentile")
        p.toolbar.logo = None
        p.quad(left=0, right=threshold, bottom=0, top=1, color="#1f77b4", level="underlay", alpha=0.2)
        p.quad(left=threshold, right=1, bottom=0, top=1, color="#ff7f0e", level="underlay", alpha=0.2)
        p.add_layout(Label(y=0.05, x=threshold + 0.01, text='High Risk', text_font_size='8pt',
                           text_color='#ff7f0e', text_alpha=0.8, text_align='left'))
        p.add_layout(Label(y=0.05, x=threshold - 0.01, text='Low Risk', text_font_size='8pt',
                           text_color='#1f77b4', text_alpha=0.8, text_align='right'))
        scatter = p.circle('Percentile', 'P(Death)', source=self.source, color="steelblue", alpha=0.1, 
                           size=7, hover_alpha=0.5, line_color=None, hover_fill_color="midnightblue",
                           hover_line_color="white", selection_color="midnightblue", selection_alpha=0.7,
                           selection_line_color="white")
        p.add_tools(HoverTool(renderers=[scatter], mode='vline', tooltips=None))
        return p

    def create_histogram_plot(self):
        p = figure(title='AML Epigenomic Risk - Select clusters on the map and their prognosis will appear here', 
                   width=self.width, height=300, x_axis_label='P(Death)', y_axis_label='Frequency', tools="save")
        p.toolbar.logo = None
        hist, edges = np.histogram(self.df2['P(Death)'], bins=50, range=[0, 1])
        hist_source = ColumnDataSource(data=dict(top=hist, left=edges[:-1], right=edges[1:]))
        p.quad(top='top', bottom=0, left='left', right='right', source=hist_source, fill_color="navy", line_color="white", alpha=0.5)
        avg_p_death = self.df2['P(Death)'].mean()
        p.add_layout(Span(location=avg_p_death, dimension='height', line_color='black', line_dash='dashed', line_alpha=0.8))
        p.add_layout(Label(x=avg_p_death, y=max(hist), text=f'Avg: {avg_p_death:.2f}', text_font_size='8pt', text_color='black', text_align='center'))
        return p, hist_source, edges

    def create_race_plot(self):
        race_counts = self.df2['Race or ethnic group'].value_counts(normalize=True) * 100
        race_hist_source = ColumnDataSource(data=dict(race=list(race_counts.index), counts=race_counts))
        p = figure(title='Race or Ethnic Group', width=self.width, height=300,
                   x_range=FactorRange(*race_counts.index), y_axis_label='Percentage', y_range=(0, 100), tools="save")
        p.toolbar.logo = None
        p.add_tools(HoverTool(tooltips=[('Count', '@counts{0.0}%')]))
        p.vbar(x='race', top='counts', width=0.9, source=race_hist_source, fill_color="green", line_color="white", alpha=0.5)
        return p, race_hist_source

    def update_histogram(self, hist_source, edges, p3):
        callback = CustomJS(args=dict(source=self.source, hist_source=hist_source, edges=edges, p3=p3), code="""
            const indices = source.selected.indices;
            const data = source.data;
            const hist_data = hist_source.data;
            const hist = new Array(edges.length - 1).fill(0);
            let sum_p_death = 0;
            for (let i = 0; i < indices.length; i++) {
                const value = data['P(Death)'][indices[i]];
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
            const avg_line = p3.renderers.find(r => r.name === 'avg_line');
            avg_line.location = avg_p_death;
            const avg_label = p3.renderers.find(r => r.name === 'avg_label');
            avg_label.x = avg_p_death;
            avg_label.text = `Avg: ${avg_p_death.toFixed(2)}`;
            p3.request_render();
        """)
        self.source.selected.js_on_change('indices', callback)

    def update_race_histogram(self, race_hist_source):
        callback = CustomJS(args=dict(source=self.source, race_hist_source=race_hist_source), code="""
            const indices = source.selected.indices;
            const data = source.data;
            const race_hist_data = race_hist_source.data;
            const race_counts = {};
            for (let race of race_hist_data['race']) {
                race_counts[race] = 0;
            }
            for (let i = 0; i < indices.length; i++) {
                const race = data['Race or ethnic group'][indices[i]];
                if (race in race_counts) {
                    race_counts[race] += 1;
                }
            }
            const total_selected = indices.length;
            for (let i = 0; i < race_hist_data['race'].length; i++) {
                race_hist_data['counts'][i] = (total_selected > 0) ? (race_counts[race_hist_data['race'][i]] / total_selected) * 100 : 0;
            }
            race_hist_source.change.emit();
        """)
        self.source.selected.js_on_change('indices', callback)

    def create_scatter_plot(self, col, test_sample=None):
        factors = [str(val) for val in self.df2[col].unique() if pd.notnull(val)]
        color_mapper = CategoricalColorMapper(factors=factors, palette=self.custom_color_palette)
        p = figure(title='Acute Leukemia Methylome Atlas', width=self.width, height=600,
                   tools="pan,wheel_zoom,box_select,reset,save", tooltips=[(str(col), '@{' + str(col) + '}')], 
                   x_axis_label=self.xaxis, y_axis_label=self.yaxis, active_drag="box_select", 
                   x_range=self.x_range, y_range=self.y_range)
        p.toolbar.logo = None
        p.toolbar_location = 'above'
        p.xaxis.axis_label_text_font_size = "8pt"
        p.yaxis.axis_label_text_font_size = "8pt"
        p.xaxis.axis_label_text_font_style = "normal"
        p.yaxis.axis_label_text_font_style = "normal"
        for factor in factors:
            view = CDSView(filter=GroupFilter(column_name=col, group=factor))
            p.scatter(x=self.xaxis, y=self.yaxis, source=self.source, view=view, 
                      color={'field': col, 'transform': color_mapper}, size=3, alpha=0.8, radius=0.2)
        if test_sample:
            test_x, test_y = self.df2.loc[test_sample][[self.xaxis, self.yaxis]]
            p.renderers.extend([
                Span(location=test_x, dimension='height', line_color="black", line_dash='dashed', line_alpha=0.8),
                Span(location=test_y, dimension='width', line_color="black", line_dash='dashed', line_alpha=0.8)
            ])
            p.star(x=test_x, y=test_y, size=15, color="black", alpha=0.9,
                   legend_label=f'Sample: {test_sample}\nPrediction: {self.df2.loc[test_sample]["AL Epigenomic Subtype"]}',
                   line_color="black", line_width=1)
            p.legend.click_policy = "hide"
        legend_items = [LegendItem(label=factor, renderers=[r]) for factor, r in zip(factors, p.renderers)]
        p.add_layout(Legend(items=legend_items, location="top", click_policy="hide", label_text_font_size="8pt", 
                            label_text_font_style="normal", glyph_height=15, glyph_width=15, spacing=1), 'right')
        return TabPanel(child=p, title=col)

    def plot(self, table=True, test_sample=None):
        data_table = self.create_data_table()
        risk_plot = self.create_risk_plot()
        hist_plot, hist_source, edges = self.create_histogram_plot()
        race_plot, race_hist_source = self.create_race_plot()

        if test_sample:
            self.update_test_sample(risk_plot, hist_plot, test_sample)

        self.update_histogram(hist_source, edges, hist_plot)
        self.update_race_histogram(race_hist_source)

        tabs = [self.create_scatter_plot(col, test_sample) for col in self.cols]
        layout = column(Tabs(tabs=tabs, tabs_location='above'), hist_plot, risk_plot, data_table) if table else column(Tabs(tabs=tabs, tabs_location='above'), hist_plot, risk_plot)
        show(layout)

    def update_test_sample(self, risk_plot, hist_plot, test_sample):
        test_val = self.df2.loc[test_sample]['P(Death)']
        vline = Span(location=test_val, dimension='height', line_color='black', line_dash='dashed', line_alpha=0.8)
        risk_plot.renderers.extend([vline])
        hist_plot.renderers.extend([vline])
        risk_plot.star(x=test_val, y=0, size=15, color="black", alpha=0.9, 
                       legend_label=f'{test_sample}, {self.df2.loc[test_sample]["AML Epigenomic Risk"]} Epigenomic Risk ({test_val:.2f})',
                       line_color="black", line_width=1)
        hist_plot.star(x=test_val, y=max(hist_plot.renderers[0].data_source.data['top']) * 0.5, size=15, color="black", alpha=0.9, 
                       legend_label=f'{test_sample}, {self.df2.loc[test_sample]["AML Epigenomic Risk"]} Epigenomic Risk ({test_val:.2f})',
                       line_color="black", line_width=1)
        risk_plot.legend.location = "bottom_right"
        risk_plot.legend.click_policy = "hide"

# Usage
# plotter = Linked_ALMA_Plotter(df)
# plotter.plot(table=True, test_sample=None)
