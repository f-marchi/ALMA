from bokeh.models import (Slider, TabPanel, Tabs, Legend, ColumnDataSource, LegendItem,
                          CDSView, GroupFilter, CategoricalColorMapper, Label, Span,
                          DataTable, HoverTool, TableColumn)
from bokeh.io import curdoc, show
from bokeh.plotting import figure
from bokeh.embed import json_item
from bokeh.layouts import column
from bokeh.themes import Theme
import pandas as pd
import json

import numpy as np
import pacmap
import pickle
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

cols = ['PaCMAP Output','Hematopoietic Entity','WHO 2022 Diagnosis','ELN 2022 Diagnosis', 'FAB', 'FLT3 ITD', 'Age (group years)',
        'Complex Karyotype', 'Primary Cytogenetic Code' , 'Sex', 'MRD 1 Status',
        'Leucocyte counts (10⁹/L)', 'Risk Group', 'Race or ethnic group',
        'Clinical Trial', 'Vital Status','First Event','Sample Type', 'Train Test']

def setup_source_dataset(pacmap_path):
    # Load clinical data
    discovery_clinical_data = pd.read_csv('../public/discovery_clinical_data.csv',
                                        low_memory=False, index_col=0)

    # Load clinical data
    validation_clinical_data = pd.read_csv('../public/validation_clinical_data.csv',
                                            low_memory=False, index_col=0)

    # Adjust clinical data
    discovery_clinical_data['Train Test'] = 'Discovery (train) Samples'
    validation_clinical_data['Train Test'] = 'Validation (test) Samples'

    discovery_clinical_data['PaCMAP Output'] = 'Patient Samples'
    validation_clinical_data['PaCMAP Output'] = 'Patient Samples'

    # Set the theme for the plot
    curdoc().theme = 'light_minimal' # or 'dark_minimal'

    df2 = pd.read_csv(pacmap_path, index_col=0)

    # Concatenate discovery and validation clinical data
    clinical_data = pd.concat([discovery_clinical_data, validation_clinical_data])

    # Join clinical data to the embedding
    df2 = df2.join(clinical_data[cols], rsuffix='_copy', on='index')

    methylscore_df = pd.read_excel('../public/ewas_cog_os_MethylScoreAML_Px.xlsx', index_col=0)

    # Concatenate df2 with temp based on the index
    df3 = df2.join(methylscore_df[['WHO 2022 Diagnosis', 'Vital Status', 'First Event', 'Clinical Trial', 'MethylScoreAML_Px', 'FLT3 ITD', 'Gene Fusion']], rsuffix='_copy', on='index').dropna(subset=['MethylScoreAML_Px'])

    df3 = df3.sort_values(by=['MethylScoreAML_Px'])

    df3 = df3.reset_index(drop=True)  

    # Normalize the Patient Number between 0 and 100 in a new column called "Percentile"
    df3['Percentile'] = df3.index / (len(df3.index) - 1)

    # Concatenate MethylScore to 2 decimal places
    df3['MethylScoreAML_Px'] = df3['MethylScoreAML_Px'].round(2)

    return df3

def create_risk_plot(df3, export=True):
    # Create a ColumnDataSource from df: source
    source = ColumnDataSource(df3)

    p1 = figure(width=1175, height=450, 
                tools="pan,wheel_zoom,box_zoom,xbox_select,reset", 
                active_drag="xbox_select")
    p1.toolbar.logo = None

    BINARY_THRESHOLD = 0.2208

    # Set the x_range and y_range of the plot to (-3, 3) and (-0.05, 1.1) respectively
    p1.x_range.bounds = (-3, 3)
    p1.x_range.start = -3  
    p1.x_range.end = 3
    p1.y_range.bounds = (-0.05, 1.1)
    p1.y_range.start = -0.05
    p1.y_range.end = 1.1

    # Add circle glyphs to the figure p with the selected and non-selected properties
    p1.circle(x='MethylScoreAML_Px', y='Percentile', source=source, selection_color='#ff7f0e',
              nonselection_alpha=1.0, color='#1f77b4', size=5, alpha=0.8, hover_color='#ff7f0e', hover_alpha=1.0)

    # Axis labels and range
    p1.xaxis.axis_label = "MethylScore"
    p1.yaxis.axis_label = "est. probability of event"
    p1.xaxis.axis_label_text_font_size = "10pt"
    p1.yaxis.axis_label_text_font_size = "10pt"
    p1.xaxis.axis_label_text_font_style = "normal"
    p1.yaxis.axis_label_text_font_style = "normal"

    # Add vertical line at the binary threshold
    p1.add_layout(Span(location=BINARY_THRESHOLD, dimension='height', line_color='black'))

    # Text annotations for risk categories
    labels = [
        Label(x=-BINARY_THRESHOLD, y=20, text="Low Risk", text_color='black', text_align='center', y_units='screen'),
        Label(x=BINARY_THRESHOLD + 0.45, y=20, text="High Risk", text_color='black', text_align='center', y_units='screen'),
    ]

    for label in labels:
        p1.add_layout(label)

    # Existing columns
    columns = [
        TableColumn(field="WHO 2022 Diagnosis", title="WHO 2022 Diagnosis"),
        TableColumn(field="Vital Status", title="Vital Status"),
        TableColumn(field="First Event", title="First Event"),
        TableColumn(field="Clinical Trial", title="Clinical Trial"),
        TableColumn(field="MethylScoreAML_Px", title="MethylScoreAML_Px"),
        TableColumn(field="FLT3 ITD", title="FLT3 ITD"),
        TableColumn(field="Gene Fusion", title="Gene Fusion")
    ]

    data_table = DataTable(source=source, columns=columns, width=1175, editable=True)


    # Define the tooltips for the HoverTool
    tooltips = [
        ("MethylScore", "@MethylScoreAML_Px"),
        ("Percentile", "@Percentile"),
    ]
    hover = HoverTool(tooltips=tooltips)
    p1.add_tools(hover)

    #### ------------------ PaCMAP Plot ------------------ ####

    # Custom color palette
    custom_color_palette = [
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
        '#a2769e',  # Soft purple
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

    # Initial setup
    title = ''
    x_range = (-45, 45)
    y_range = (-45, 45)
    datapoint_size = 5
    tooltip_dx_cols = 'WHO 2022 Diagnosis'
    width = 1000
    height = 1000

    # Initialize Bokeh Document with theme
    curdoc().theme = Theme(json={
        "attrs": {
            "Axis": {
                "major_label_text_color": 'black',
                "major_label_text_font": 'Arial',
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

    tabs = []
    slider = Slider(title="Adjust datapoint size", start=0, end=10, step=1, value=datapoint_size)

    for col in cols:
        factors = [str(val) for val in df3[col].unique() if pd.notnull(val)]
        color_mapper = CategoricalColorMapper(factors=factors, palette=custom_color_palette)

        p = figure(title=title, width=width, height=height, x_range=x_range, y_range=y_range, 
                tools="pan,wheel_zoom,reset,save,box_select", tooltips=[("Dx", "@{"+tooltip_dx_cols+"}")], 
                x_axis_label='Longitude (PaCMAP 1)', y_axis_label='Latitude (PaCMAP 2)', active_drag="box_select")
        p.toolbar.logo = None
        p.toolbar_location = 'above'
        p.xaxis.axis_label_text_font_size = "10pt"
        p.yaxis.axis_label_text_font_size = "10pt"
        p.xaxis.axis_label_text_font_style = "normal"
        p.yaxis.axis_label_text_font_style = "normal"

        # Create scatter plot for each factor
        for factor in factors:
            view = CDSView(filter=GroupFilter(column_name=col, group=factor))
            p.scatter(x="PaCMAP 1", y="PaCMAP 2", source=source, view=view, 
                    color={'field': col, 'transform': color_mapper})

        # Create a list of legend items
        legend_items = [LegendItem(label=factor, renderers=[r]) for factor, r in zip(factors, p.renderers)]

        # Create a legend
        legend = Legend(items=legend_items, location="center")

        # Add the legend to the plot
        p.add_layout(legend, 'below')

        tab = TabPanel(child=p, title=col)
        tabs.append(tab)

    tabs_control = Tabs(tabs=tabs, tabs_location='left')
    layout = column(p1, data_table, tabs_control, slider)

    if (export):
        # Export as BokehJS
        json_layout = json.dumps(json_item(layout, "myplot"))
        with open("bokeh_output.json", 'w') as f:
            f.write(json_layout)

    return layout

def bed_to_pacmap(file_path, plot_filename):
    # Column names for the two parts
    columns_part1 = ["chrom", "start_position", "end_position", "modified_base_code", "score", "strand", "start_position_compat", "end_position_compat", "color"]
    columns_part2 = ["N_valid_cov", "fraction_modified", "N_mod", "N_canonical", "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall"]

    print(f'Processing nanopore...')

    # Read the first part with '\t' delimiter
    df_part1 = pd.read_csv(file_path, sep='\t', header=None, usecols=range(9), names=columns_part1)

    # Read the second part with ' ' delimiter
    df_part2 = pd.read_csv(file_path, delim_whitespace=True, header=None, usecols=range(9, 18), names=columns_part2)

    # Concatenate the two dataframes horizontally
    df = pd.concat([df_part1, df_part2], axis=1)

    # Filter the dataframe by modified base code
    df_5mc = df[df["modified_base_code"] == "m"]

    # Copy the df_5mc dataframe to avoid SettingWithCopyWarning
    df_5mc_copy = df_5mc.copy()

    print(f'Processing list of suboptimal probes...')
    zhou2016_probes = pd.read_csv('../public/EPIC.anno.GRCh38.tsv', sep='\t')

    # Add a column with the coordinate of the probe in the copied dataframe
    df_5mc_copy['coordinate'] = df_5mc_copy['chrom'].astype(str) + ':' + df_5mc_copy['start_position'].astype(str)
    zhou2016_probes['coordinate'] = zhou2016_probes['chrm'].astype(str) + ':' + zhou2016_probes['start'].astype(str)

    # Merge the two dataframes
    df_5mc_merged = pd.merge(df_5mc_copy, zhou2016_probes[['probeID','coordinate']], on='coordinate', how='inner')

    # Set the index to the probeID
    df_5mc_merged = df_5mc_merged.set_index('probeID')

    # Create beta values column
    df_5mc_merged['Kasumi1-naive'] = (df_5mc_merged['fraction_modified'] / 100).round(3)

    # Create a new dataframe with only the beta values
    df_nanopore = df_5mc_merged[['Kasumi1-naive']].T

    # # Create a new dataframe to contain metadata for `df_nanopore`, which represents sample 'UF_hem_1832_PB'
    # validation_clinical_data = pd.DataFrame(index=df_nanopore.index)
    # discovery_clinical_data = pd.read_csv('../public/discovery_clinical_data.csv', index_col=0)

    # # Adjust clinical data
    # discovery_clinical_data['Train Test'] = 'Discovery (train) Samples'
    # validation_clinical_data['Train Test'] = 'Kasumi1-naive'

    # discovery_clinical_data['PaCMAP Output'] = 'Patient Samples'
    # validation_clinical_data['PaCMAP Output'] = 'Patient Samples'

    # discovery_clinical_data['Batch'] = df_discovery['Batch']
    # validation_clinical_data['Batch'] = 'Nanopore'

    # # add the following columns to `validation_clinical_data`: 'Clinical Trial', 'Sample Type', 'Patient_ID', 'ELN AML 2022 Diagnosis', 'Hematopoietic Lineage'
    # validation_clinical_data['Clinical Trial'] = 'UF Hem Bank'
    # validation_clinical_data['Sample Type'] = 'Peripheral Blood'
    # validation_clinical_data['Patient_ID'] = 'Kasumi1-naive'
    # validation_clinical_data['ELN AML 2022 Diagnosis'] = np.nan
    # validation_clinical_data['Hematopoietic Lineage'] = np.nan

    print(f'Processing discovery CpGs...')
    df_discovery = pd.read_pickle('../public/3308samples_333059cpgs_withbatchcorrection_bvalues.pkl.gz').sort_index()
    
    print(f'Fetching common CpGs between discovery and sample...')

    # use overlapping features between df_discovery and df_nanopore
    common_features = [x for x in df_discovery.columns if x in df_nanopore.columns]

    # apply `common_features` to both df_discovery and df_nanopore
    df_discovery = df_discovery[common_features]
    df_nanopore = df_nanopore[common_features]

    print(
    f' Discovery dataset (df_discovery) contains {df_discovery.shape[1]} \
    columns (5mC nucleotides/probes) and {df_discovery.shape[0]} rows (samples).')

    print(
    f' Validation dataset (df_nanopore) contains {df_nanopore.shape[1]} \
    columns (5mC nucleotides/probes) and {df_nanopore.shape[0]} rows (samples).')

    print(f'Running PaCMAP...')

    reducer = pacmap.PaCMAP(n_components=2, n_neighbors=15, MN_ratio=0.4, FP_ratio=16.0, 
                        random_state=42, lr=0.1, num_iters=5000)

    embedding_training = reducer.fit_transform(df_discovery.to_numpy(dtype='float16'))
    embedding_validation = reducer.transform(df_nanopore.to_numpy(dtype='float16'), df_discovery.to_numpy(dtype='float16'))

    # Plot embedding
    plt.figure(figsize=(10, 10))
    plt.scatter(embedding_training[:, 0], embedding_training[:, 1], c='green', s=5, alpha=0.5)
    plt.scatter(embedding_validation[:, 0], embedding_validation[:, 1], c='black', s=50, alpha=1)
    plt.xticks([])
    plt.yticks([])

    print(f'Done!')

    # Save the plot
    plt.savefig(plot_filename, dpi=300)