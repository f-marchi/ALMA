"""
This module implements custom functions for a combined Bokeh plot with two scatter plots and a table linked.

"""  
import numpy as np
import pandas as pd
from bokeh.layouts import column, gridplot
from bokeh.plotting import figure, show, save, output_notebook
from bokeh.models import (TabPanel, Tabs, Legend, ColumnDataSource, LegendItem,
                          CDSView, GroupFilter, CategoricalColorMapper, Label,
                          DataTable, HoverTool, TableColumn, Span)
from sklearn.metrics import roc_curve, auc



def get_custom_color_palette():
    list = [
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
    return list

def plot_linked_scatters(df, table=True, test_sample=None, 
                        xaxis = "PaCMAP 1 of 2", yaxis = "PaCMAP 2 of 2",
                        x_range= (-45,40), y_range=(-50,45), 
                        cols = ['AL Epigenomic Subtype', 'Hematopoietic Entity','WHO 2022 Diagnosis', 
                                  'Vital Status','AML Epigenomic Risk', 'Risk Group AAML1831', 'Clinical Trial',
                                  'Race or ethnic group','Age (group years)'],
                        p_death = 'P(Death) at 5y'):

    # Rank samples by P(Death) and call it "Percentile"
    df_px = df[~df[p_death].isna()]
    df_px2 = df_px.sort_values(by=p_death).reset_index().reset_index(names=['Percentile']).set_index('index')
    df_px2['Percentile'] = df_px2['Percentile'] / len(df_px2['Percentile'])
    df2 = df.join(df_px2[['Percentile']])
    
    source = ColumnDataSource(df2)
    width = 1000
    font_size = "8pt"
    x= p_death
    y = 'Percentile'
    threshold = 0.5

    custom_color_palette = get_custom_color_palette()

    columns = [TableColumn(field=col, title=col) for col in cols]
    columns.append(TableColumn(field='Gene Fusion', title='Gene Fusion'))
    columns.append(TableColumn(field='Karyotype', title='Karyotype'))


    data_table = DataTable(source=source, columns=columns, editable=True, width=width,
                        index_position=None, height=300)

    p1 = figure(title= 'AML Epigenomic Risk' ,width=width, height=300,
                tools="xbox_select,reset,save", active_drag='xbox_select',
                x_axis_label=x, y_axis_label="est. probability of progression")
    p1.toolbar.logo = None

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

    scatter1 = p1.circle(y, x, source=source, color="steelblue", alpha=0.1, 
                size=7, hover_alpha=0.5, line_color=None, hover_fill_color="midnightblue",
                hover_line_color="white", selection_color="midnightblue", selection_alpha=0.7,
                selection_line_color="white", )

    scatter1_hover_tool = HoverTool(renderers=[scatter1], mode='vline', tooltips=None)

    if test_sample:
        vline = Span(location=df2.loc[test_sample][x],
                     dimension='height', line_color='black', line_dash='dashed', line_alpha=0.8)
        p1.renderers.extend([vline])

        p1.star(x=df2.loc[test_sample][x],
                y=df2.loc[test_sample][y],
                size=15, color="black", alpha=0.9, 
                legend_label=f'{test_sample}, {df2.loc[test_sample]["AML Epigenomic Risk"]} Epigenomic Risk ({df2.loc[test_sample][x]:.2f})',
                line_color="black", line_width=1)

        # move `p1.star` legend to bottom right
        p1.legend.location = "bottom_right"
        p1.legend.click_policy = "hide"

    tabs = []

    for col in cols:
        factors = [str(val) for val in df2[col].unique() if pd.notnull(val)]
        color_mapper = CategoricalColorMapper(factors=factors, palette=custom_color_palette)

        p2 = figure(title='Acute Leukemia Methylome Atlas', width=width, height=600,
                    tools="pan,wheel_zoom,box_select,reset,save", tooltips= [(str(col),'@{'+ str(col)+'}')], 
                    x_axis_label=xaxis, y_axis_label=yaxis,
                    active_drag="box_select", x_range= x_range, y_range=y_range)

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
            p2.scatter(x=xaxis, y=yaxis, source=source, view=view, 
                    color={'field': col, 'transform': color_mapper}, size=3, alpha=0.8, radius=0.2)
        if test_sample:
                vline = Span(location=df2.loc[test_sample][xaxis],
                            dimension='height', line_color="black", line_dash='dashed',line_alpha=0.8)
                hline = Span(location=df2.loc[test_sample][yaxis],
                            dimension='width', line_color="black", line_dash='dashed',line_alpha=0.8)
                p2.renderers.extend([vline, hline])
                p2.star(x=df2.loc[test_sample][xaxis], y=df2.loc[test_sample][yaxis],
                size=15, color="black", alpha=0.9, legend_label=f'Sample: {test_sample}\nPrediction: {df2.loc[test_sample]["AL Epigenomic Subtype"]}',
                 line_color="black", line_width=1)
                p2.legend.click_policy = "hide"

        # Create a list of legend items
        legend_items = [LegendItem(label=factor, renderers=[r]) for factor, r in zip(factors, p2.renderers)]

        # Create a legend
        legend = Legend(items=legend_items, location="top", click_policy="hide",
                        label_text_font_size=font_size, label_text_font_style="normal",
                        glyph_height=15, glyph_width=15, spacing=1, 
                        # margin=5, padding=0,
                        )

        # Add the legend to the plot
        p2.add_layout(legend, 'right')

        tab = TabPanel(child=p2, title=col)
        tabs.append(tab)

    tabs_control = Tabs(tabs=tabs, tabs_location='above')

    p1.add_tools(scatter1_hover_tool)

    # save(column(tabs_control,p1,data_table), filename="FM_interactive_plots_1-30-24_preliminary.html", title="AL Methylome Atlas")
    if table:
        return(show(column(tabs_control,p1,data_table)))
    else:
        return(show(column(tabs_control,p1)))

def plot_roc_auc_with_riskgroup(df, target, model_name, risk_group='Risk Group', title=None, sum_models=False):
    """
    Plots ROC AUC flexibly using Bokeh.
    
    Parameters:
    - df: pandas DataFrame containing model predictions as columns and actual target variable.
    - target: Name of the column containing the actual target variable.
    - model_name: Name of the column containing the model predictions.
    - risk_group: Name of the column containing the risk group.
    - title: Title of the plot.
    """

    def category_to_integer(df, model_name, risk_group=None, sum_models=sum_models):

        df_ = df.copy()        
        low_high_dict = {'Low': 0, 'Low Risk': 0,
                        'Standard':0.5, 'Standard Risk': 0.5,
                        'High': 1, 'High Risk': 1}

        if df[model_name].dtype == 'O':
            df_[model_name] = df_[model_name].map(low_high_dict)
            
        df_[risk_group] = df_[risk_group].map(low_high_dict)

        if sum_models:
            df_[model_name + ' + ' + risk_group] = (df_[model_name] + df_[risk_group])/2
            df_ = df_[[model_name + ' + ' + risk_group, target]]
        else:
            df_ = df_[[model_name, risk_group, target]]

        # drop rows with missing values
        df_ = df_.dropna()

        return df_

    df = category_to_integer(df, model_name, risk_group=risk_group)

    # colors = itertools.cycle(Spectral11)
    colors = ['navy', 'firebrick', 'olive']

    if title:
        title_ = title + ', n=' + str(len(df))
    else:
        title_ = ''

    p = figure(title=title_,
               x_axis_label='False Positive Rate',
               y_axis_label='True Positive Rate',
               width=325, height=325,
               tools='save,reset,pan')
    
    p.line([0, 1], [0, 1], line_dash="dashed", color="gray", line_width=1)

    for column, color in zip(df.columns.difference([target]), colors):
        fpr, tpr, _ = roc_curve(df[target], df[column])
        roc_auc = auc(fpr, tpr)
        p.line(fpr, tpr, legend_label=f"{column}\nAUC = {roc_auc:.2f}",
               color=color, line_width=2, alpha=0.8)

    p.legend.location = "bottom_right"
    p.legend.click_policy="hide"
    p.toolbar.logo = None
    p.legend.label_text_font_size = '8pt'
    p.legend.spacing = 2
    p.xaxis.axis_label_text_font_style = "normal"
    p.yaxis.axis_label_text_font_style = "normal"
    p.legend.background_fill_alpha = 0.8
    p.title.text_font_size = '9pt'

    return p

def plot_multiclass_roc_auc(df, target_columns, title=None):
    """
    Plots ROC AUC for multiple classes using Bokeh, targeting multiple classes,
    with legend below the plot.
    """
    
    colors = get_custom_color_palette()
    
    if title:
        title_ = title + ', n=' + str(len(df))
    else:
        title_ = ''

    # Initialize figure without a legend inside it
    p = figure(title=title_,
               x_axis_label='False Positive Rate',
               y_axis_label='True Positive Rate',
               width=350, height=800,
               tools='save,reset,pan')
    
    p.line([0, 1], [0, 1], line_dash="dashed", color="gray", line_width=1)

    legend_items = []  # To hold the legend items

    for target, color in zip(target_columns, colors):
        column = 'P(' + target + ')'
        if column in df.columns:  # Check if the probability column exists
            fpr, tpr, _ = roc_curve(df[target], df[column])
            roc_auc = auc(fpr, tpr)
            line = p.line(fpr, tpr,
                          color=color, line_width=2, alpha=0.8)
            # Create legend items manually
            legend_items.append(LegendItem(label=f"{column} ({roc_auc:.2f})", renderers=[line]))

    # Create a Legend with the collected items and add it below the plot
    legend = Legend(items=legend_items, spacing=0, location="left",
                    label_text_font_size='7pt', click_policy="hide",
                    background_fill_alpha=0.8)
    p.add_layout(legend, 'below')


    p.toolbar.logo = None
    p.xaxis.axis_label_text_font_style = "normal"
    p.yaxis.axis_label_text_font_style = "normal"
    p.xaxis.axis_label_text_font_size = '8pt'
    p.yaxis.axis_label_text_font_size = '8pt'
    p.title.text_font_size = '10pt'

    return p

def process_dataset_for_multiclass_auc(df):
    # One hot encode `df_dx['AL Epigenomic Subtype']`
    df_dx_dummies = pd.get_dummies(df['WHO 2022 Diagnosis'])

    # transform boolean columns to integer
    df_dx_dummies = df_dx_dummies.astype(int)

    # join the one hot encoded columns with the original dataframe
    df_dx_auc = pd.concat([df.iloc[:, -26:-1], df_dx_dummies], axis=1)

    return df_dx_auc, df_dx_dummies