#!/usr/bin/env python
# coding: utf-8

# # DNA methylation-based diagnosis and prognosis of pediatric AML

# We propose to leverage machine learning tools to develop DNA methylation-based signatures of clinical utility in pediatric AML.

# In[1]:


import pandas as pd

PaCMAP_path = '../Data/Processed_Data/PaCMAP_Results/'
input_path = '../Data/Processed_Data/'
output_path = '../Data/Processed_Data/'

x_train = pd.read_pickle(PaCMAP_path+'embedding.pkl')
x_test = pd.read_pickle(PaCMAP_path+'embedding_test.pkl')

y = pd.read_csv(input_path+'y.csv', index_col=0)

labels = pd.read_excel(input_path+'y_plus_WHOclass.xlsx', index_col=0)
labels = pd.concat([y, labels], axis=1)
labels = labels[labels.index.isin(x_train.index)]['WHO Classification']

y = y.join(labels.to_frame('WHO Classification'))
y['KMT2A Fusions'] = y[y['WHO Classification'].isin(['AML with KMT2A-rearrangement'])]['Gene Fusion'] 

y_train = y[~y['Clinical Trial'].isin(['AML02','AML08'])]
y_test = y[y['Clinical Trial'].isin(['AML02','AML08'])]

from bokeh.plotting import figure, show
from bokeh.transform import factor_cmap

from bokeh.layouts import layout
from bokeh.models import Div, Slider, TabPanel, Tabs
from bokeh.io import curdoc, output_notebook

output_notebook()


# In[2]:


list = ['Primary Cytogenetic Code', 'FAB', 'FLT3 ITD','Age group (years)',
       'WHO Classification','Complex Karyotype', 'Karyotype']

df = x_train.join(y_train[list]).reset_index() # join embedding with labels
df['PaCMAP Output'] = 'PaCMAP Output'

curdoc().theme = 'light_minimal'

def fig():
    """
    Figure specs for Bokeh plot
    """
    
    fig = figure(
           width=600,
           height=600,
           sizing_mode='fixed',
           x_axis_label='PaCMAP 1',
           y_axis_label='PaCMAP 2',
           tools="pan,wheel_zoom, reset, save",
           active_drag="pan",
           active_scroll="wheel_zoom",
           tooltips=[("Sample", "@index"),
                     ("Karyotype", "@Karyotype"),])

    return(fig)

def scatter(df, p, hue):
    """
    Scatter plot of embedding with color by hue
    
    Parameters
    ----------
    p : bokeh.plotting.figure.Figure
        Bokeh figure object
    hue : str
        Column name of df to color by
    Returns
    -------
    points : bokeh.models.renderers.GlyphRenderer
        Bokeh glyph renderer object
        
    """
    df = df[~df[hue].isna()] # df where df hue is not nan
    points = p.scatter(x="PaCMAP 1",
                   y= "PaCMAP 2",
                   source=df.copy(),
                   fill_alpha=0.8,
                   size=5,
                   color=factor_cmap(field_name= hue,
                                     palette='Category10_10',
                                     factors= sorted(df[hue].unique())),
                    legend_group=hue)
    return(points)

p1 = fig()
points1 = scatter(df, p1, hue='PaCMAP Output')
tab1 = TabPanel(child=p1, title='PaCMAP Output')
p1.toolbar.logo = None

p2 = fig()
points2 = scatter(df, p2, hue='FAB')
tab2 = TabPanel(child=p2, title="FAB")
p2.toolbar.logo = None

p3 = fig()
points3 = scatter(df, p3, hue='Complex Karyotype')
tab3 = TabPanel(child=p3, title="Complex Karyotype")
p3.toolbar.logo = None

p4 = fig()
points4 = scatter(df, p4, hue='FLT3 ITD')
tab4 = TabPanel(child=p4, title="FLT3 ITD")
p4.toolbar.logo = None

p5 = fig()
points5 = scatter(df, p5, hue='Primary Cytogenetic Code')
tab5 = TabPanel(child=p5, title='Primary Cytogenetic Code')
p5.toolbar.logo = None

p6 = fig()
points6 = scatter(df, p6, hue='WHO Classification')
tab6 = TabPanel(child=p6, title='WHO Classification')
p6.toolbar.logo = None
# Remove legend
p6.legend.visible = False

p7 = fig()
points7 = scatter(df, p7, hue='Age group (years)')
tab7 = TabPanel(child=p7, title='Age group (years)')
p7.toolbar.logo = None


tabs = Tabs(tabs=[tab1, tab2, tab3, tab4, tab5, tab6, tab7], tabs_location='left')

div = Div(
    text="""
          <b> The AML Diagnostic Map</b>
            <br> Interactive visualization of the pediatric AML methylome:</br>
          """,
    width=200,
    height=85)

slider = Slider(
    title="Adjust datapoint size",
    start=0,
    end=20,
    step=1,
    value=(points1.glyph.size))

slider.js_link("value", points1.glyph, "size")
slider.js_link("value", points2.glyph, "size")
slider.js_link("value", points3.glyph, "size")
slider.js_link("value", points4.glyph, "size")
slider.js_link("value", points5.glyph, "size")
slider.js_link("value", points6.glyph, "size")

# create layout
layout = layout([[[div,tabs, slider]]])

# show result
show(layout)


# ## Table of Contents
# 
# ```{tableofcontents}
# ```
# 
