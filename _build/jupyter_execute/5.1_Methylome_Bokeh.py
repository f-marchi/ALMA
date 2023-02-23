#!/usr/bin/env python
# coding: utf-8

# # PaCMAP Interactive Visualization

# ## Where the data at?

# In[13]:


input_path = '../Data/Processed_Data/'
output_path = '../Data/Processed_Data/'


# ## Load Datasets

# In[14]:


import pandas as pd

x_train = pd.read_pickle(input_path+'embedding.pkl')
x_test = pd.read_pickle(input_path+'embedding_test.pkl')

nanopore_sample = pd.read_pickle(input_path+'embedding_nano.pkl')

y = pd.read_csv(input_path+'y.csv', index_col=0)

labels = pd.read_excel('../Data/Raw_Data/Clinical_Data/labelsCOG_WHOClass.xlsx', index_col=0)
labels = pd.concat([y, labels], axis=1)
labels = labels[labels.index.isin(x_train.index)]['WHO Classification']
y = y.join(labels.to_frame('WHO Classification'))

y_train = y[~y['Clinical Trial'].isin(['AML02','AML08'])]
y_test = y[y['Clinical Trial'].isin(['AML02','AML08'])]


# ## Bokeh

# ### Import Libraries

# In[15]:


from bokeh.plotting import figure, show, output_file
from bokeh.transform import factor_cmap


# In[16]:


from bokeh.layouts import layout
from bokeh.models import Div, Slider, Spinner, ResetTool, TabPanel, Tabs
from bokeh.io import curdoc, output_notebook

output_notebook()


# ### Define Variables

# In[17]:


# p = figure(title = "AML Methylome Map", tooltips=[("Vital Status", "@{Vital Status}")])
# p.scatter("PaCMAP 1", "PaCMAP 2", source=df, fill_alpha=0.8, size=5,
#           color=factor_cmap(field_name= hue,
#                             palette='Category10_10',
#                             factors= sorted(df[hue].unique())))

# show(p)


# In[18]:


# select = Select(title="Clinical Annotation:", options=list)
# select.js_on_change("value", CustomJS(code="""
#     console.log('select: value=' + this.value, this.toString())
# """))

# show(select)



# In[19]:


list = ['Primary Cytogenetic Code', 'Vital Status', 
        'FLT3 ITD','Risk Group','Age (years)',
        'Gene Fusion', 'Sex', 'Batch', 'WHO Classification']

df = x_train.join(y_train[list]).reset_index() # join embedding with labels
#df = df[~df[hue].isna()] # df where df hue is not nan


# In[53]:


from bokeh.layouts import layout
from bokeh.palettes import Category10_10

output_file(filename="tempplot.html",
            title="AML Methylome Map")
curdoc().theme = 'light_minimal'

def fig():
    """
    Figure specs for Bokeh plot
    """
    
    fig = figure(title="A longer title\nwith a second line underneath",
           width=600,
           height=600,
           sizing_mode='fixed',
           x_axis_label='PaCMAP 1',
           y_axis_label='PaCMAP 2',
           tools="pan,wheel_zoom, reset, save",
           active_drag="pan",
           active_scroll="wheel_zoom",
           tooltips=[("Fusion Partner","@{Gene Fusion}"),
                     ("WHO Classification", "@{WHO Classification}")])

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

def draw_tab(hue):
    p1 = fig() # create figure
    points1 = scatter(df, p1, hue=hue) # create scatter plot
    tab1 = TabPanel(child=p1, title=hue) # create tab
    p1.toolbar.logo = None
    return(tab1, points1)


tab1 = draw_tab(hue='Primary Cytogenetic Code')[0]
tab2 = draw_tab(hue='Vital Status')[0]
points1 = draw_tab(hue='Primary Cytogenetic Code')[1]
points2 = draw_tab(hue='Vital Status')[1]

tabs = Tabs(tabs=[tab1, tab2])

slider = Slider(title="Adjust datapoint size", 
                    start=0,end=20,step=1,
                    value=(points1.glyph.size)) # create slider
    
slider.js_link("value", points1.glyph, "size") # link slider to glyph size
slider.js_link("value", points2.glyph, "size")

# create div
div = Div(
    text="""
          <b> The AML Diagnostic Map</b>
            <br> Interactive visualization of the pediatric AML methylome:</br>
          """,
    width=200,
    height=100,
)

# create layout
layout = layout(
    [
        [[div, slider],[tabs]]
    ]
)

# show result
show(layout)


# In[44]:


list = ['Primary Cytogenetic Code', 'Vital Status', 'Batch']

for i in range(0, len(list)):
        


# In[54]:


from bokeh.layouts import layout
#from bokeh.palettes import Category10_10

output_file(filename="tempplot.html",
            title="AML Methylome Map")
curdoc().theme = 'light_minimal'

def fig():
    """
    Figure specs for Bokeh plot
    """
    
    fig = figure(title="A longer title\nwith a second line underneath",
           width=600,
           height=600,
           sizing_mode='fixed',
           x_axis_label='PaCMAP 1',
           y_axis_label='PaCMAP 2',
           tools="pan,wheel_zoom, reset, save",
           active_drag="pan",
           active_scroll="wheel_zoom",
           tooltips=[("Fusion Partner","@{Gene Fusion}"),
                     ("WHO Classification", "@{WHO Classification}")])

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
points1 = scatter(df, p1, hue='Primary Cytogenetic Code')
tab1 = TabPanel(child=p1, title='Primary Cytogenetic Code')
p1.toolbar.logo = None

p2 = fig()
points2 = scatter(df, p2, hue='Vital Status')
tab2 = TabPanel(child=p2, title="Vital Status")
p2.toolbar.logo = None

p3 = fig()
points3 = scatter(df, p3, hue='Batch')
tab3 = TabPanel(child=p3, title="Batch")
p3.toolbar.logo = None

p4 = fig()
points4 = scatter(df, p4, hue='FLT3 ITD')
tab4 = TabPanel(child=p4, title="FLT3 ITD")
p4.toolbar.logo = None

tabs = Tabs(tabs=[tab1, tab2, tab3, tab4])

div = Div(
    text="""
          <b> The AML Diagnostic Map</b>
            <br> Interactive visualization of the pediatric AML methylome:</br>
          """,
    width=200,
    height=100)

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

# create layout
layout = layout([[[div, slider],[tabs]]])

# show result
show(layout)



# In[38]:


from bokeh.layouts import layout
#from bokeh.palettes import Category10_10

output_file(filename="tempplot.html",
            title="AML Methylome Map")
curdoc().theme = 'light_minimal'

def p():
    """
    Figure specs for Bokeh plot
    """
    
    p = figure(title="A longer title\nwith a second line underneath",
           width=600,
           height=600,
           sizing_mode='fixed',
           x_axis_label='PaCMAP 1',
           y_axis_label='PaCMAP 2',
           tools="pan,wheel_zoom, reset, save",
           active_drag="pan",
           active_scroll="wheel_zoom",
           tooltips=[("Fusion Partner","@{Gene Fusion}"),
                     ("WHO Classification", "@{WHO Classification}")])

    return(p)

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

p1 = p()
points1 = scatter(df, p1, hue='Primary Cytogenetic Code')
tab1 = TabPanel(child=p1, title='Primary Cytogenetic Code')
p1.toolbar.logo = None

p2 = p()
points2 = scatter(df, p2, hue='Vital Status')
tab2 = TabPanel(child=p2, title="Vital Status")
p2.toolbar.logo = None

p3 = p()
points3 = scatter(df, p3, hue='Batch')
tab3 = TabPanel(child=p3, title="Batch")
p3.toolbar.logo = None

p4 = p()
points4 = scatter(df, p4, hue='FLT3 ITD')
tab4 = TabPanel(child=p4, title="FLT3 ITD")
p4.toolbar.logo = None

tabs = Tabs(tabs=[tab1, tab2, tab3, tab4])

div = Div(
    text="""
          <b> The AML Diagnostic Map</b>
            <br> Interactive visualization of the pediatric AML methylome:</br>
          """,
    width=200,
    height=100)

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


# # set up Select
# select = Select(title="Clinical Annotation:",value='WHO Classification' , options=list)
# #select.js_link("value", p.color.factor_cmap, "hue")


# select.js_on_change("value", CustomJS(code="""
#      console.log('select: value=' + this.value, this.toString())
#  """))

# create layout
layout = layout([[[div, slider],[tabs]]])

# show result
show(layout)



# In[22]:


# from numpy import cos, linspace, pi, sin

# from bokeh.core.enums import LegendLocation
# from bokeh.io import show
# from bokeh.models import (Circle, ColumnDataSource, DataRange1d, Legend, Line,
#                           LinearAxis, PanTool, Plot, SaveTool, WheelZoomTool)

# x = linspace(-2*pi, 2*pi, 400)
# y = sin(x)
# y2 = cos(x)

# source = ColumnDataSource(data=dict(x=x, y=y, y2=y2))

# xdr = DataRange1d()
# ydr = DataRange1d()

# plot = Plot(
#     x_range=xdr, y_range=ydr,
#     width=1000, height=600,
#     min_border=0,
#     toolbar_location=None,
#     background_fill_color='#F0F0F0',
#     border_fill_color='lightgray',
# )

# line_glyph = Line(x="x", y="y", line_color="navy", line_width=2, line_dash="dashed")
# line = plot.add_glyph(source, line_glyph)
# circle = Circle(x="x", y="y2", size=6, line_color="red", fill_color="orange", fill_alpha=0.6)
# circle = plot.add_glyph(source, circle)

# pan = PanTool()
# wheel_zoom = WheelZoomTool()
# preview_save = SaveTool()

# plot.add_tools(pan, wheel_zoom, preview_save)

# # Add axes (Note it's important to add these before adding legends in side panels)
# plot.add_layout(LinearAxis(), 'below')
# plot.add_layout(LinearAxis(), 'left')
# plot.add_layout(LinearAxis(), 'right')

# def add_legend(location, orientation, side):
#     legend = Legend(
#         # TODO: title
#         items=[("line", [line]), ("circle", [circle])],
#         location=location, orientation=orientation,
#         border_line_color="black",
#         title='Example Title'
#     )
#     plot.add_layout(legend, side)

# # Add legends in names positions e.g. 'top_right', 'top_left' (see plot for all)
# for location in LegendLocation:
#     add_legend(location, "vertical", "center")

# # Add legend at fixed positions
# add_legend((150, 50), "horizontal", "center")

# # Add legend in side panels
# add_legend("center_left", "horizontal", "above")
# add_legend("center", "horizontal", "below")
# add_legend("center", "vertical", "left")
# add_legend("bottom_center", "vertical", "right")

# show(plot)


# In[23]:


# from bokeh.layouts import layout

# # # create plot with circle glyphs
# # p = figure(title="AML Methylome Map",
# #            width=800,
# #            height=800,
# #            sizing_mode='fixed',
# #            x_axis_label='PaCMAP 1',
# #            y_axis_label='PaCMAP 2',
# #            tools="pan,wheel_zoom, reset, save",
# #            active_drag="pan",
# #            active_scroll="wheel_zoom",
# #            tooltips=[("Vital Status", "@{Vital Status}"),
# #                      ("WHO Classification", "@{WHO Classification}")])

# # points = p.scatter(x="PaCMAP 1",
# #                    y= "PaCMAP 2",
# #                    source=df.copy(),
# #                    fill_alpha=0.8,
# #                    size=5,
# #                    color=factor_cmap(field_name= hue,
# #                                      palette='Category10_10',
# #                                      factors= sorted(df[hue].unique())))

# p = p()
# points = scatter(p, hue)

# #p.toolbar_location = None 
# p.toolbar.logo = None  
# p.title.text_font_size = '16pt'                                  

# div = Div(
#     text="""
#           <p> Explore and interact with the data:</p>
#           """,
#     width=200,
#     height=30,
# )

# # set up RangeSlider
# slider = Slider(
#     title="Adjust datapoint size",
#     start=0,
#     end=20,
#     step=1,
#     value=(points.glyph.size),
# )
# slider.js_link("value", points.glyph, "size")

# # set up Select
# select = Select(title="Clinical Annotation:",value='WHO Classification' , options=list)
# #select.js_link("value", p.color.factor_cmap, "hue")


# select.js_on_change("value", CustomJS(code="""
#      console.log('select: value=' + this.value, this.toString())
#  """))

# # create layout
# layout = layout(
#     [
#         [[p],[div, slider, select]]
#     ]
# )


# # show result
# show(layout)



# In[24]:


# from bokeh.layouts import layout
# from bokeh.models import Div, RangeSlider, Spinner
# from bokeh.plotting import figure, show

# # prepare some data
# x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# y = [4, 5, 5, 7, 2, 6, 4, 9, 1, 3]

# # create plot with circle glyphs
# p = figure(x_range=(1, 9), width=500, height=250)
# points = p.circle(x=x, y=y, size=30, fill_color="#21a7df")

# # set up textarea (div)
# div = Div(
#     text="""
#           <p>Select the circle's size using this control element:</p>
#           """,
#     width=200,
#     height=30,
# )

# # set up spinner
# spinner = Spinner(
#     title="Circle size",
#     low=0,
#     high=60,
#     step=5,
#     value=points.glyph.size,
#     width=200,
# )
# spinner.js_link("value", points.glyph, "size")

# # set up RangeSlider
# range_slider = RangeSlider(
#     title="Adjust x-axis range",
#     start=0,
#     end=10,
#     step=1,
#     value=(p.x_range.start, p.x_range.end),
# )
# range_slider.js_link("value", p.x_range, "start", attr_selector=0)
# range_slider.js_link("value", p.x_range, "end", attr_selector=1)

# # create layout
# layout = layout(
#     [
#         [div, spinner],
#         [range_slider],
#         [p],
#     ]
# )

# # show result
# show(layout)

