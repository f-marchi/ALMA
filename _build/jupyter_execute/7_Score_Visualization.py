#!/usr/bin/env python
# coding: utf-8

# # Kaplan Meiers, Forest Plots, Patient Characteristics

# ## Where the data at?

# In[1]:


input_path = '../Data/Processed_Data/'
output_path = '../Data/Processed_Data/'


# ## Load Datasets

# In[2]:


import pandas as pd

y = pd.read_csv(input_path+'y_MethylScore.csv', index_col=0)
# y = pd.read_csv(input_path+'y.csv', index_col=0)


# In[3]:


y_train = y[~y['Clinical Trial'].isin(['AML02','AML08'])]
#y_train2 = y_train[y_train['os.evnt'].notnull()]
y_test = y[y['Clinical Trial'].isin(['AML02','AML08'])]

y_test.shape, y_train.shape


# In[4]:


# Import Plotting Functions
from FM_Functions.Data_Visualization import *


# In[5]:


score_name = "MethylScore"


# ## Kaplan Meiers

# In[6]:


draw_kaplan_meier(scorename=score_name,
                        df=y_train[~y_train['Clinical Trial'].isin(['AML02','AML08','AML05'])],
                        save_plot=True,
                        add_risk_counts=False,
                        trialname='Discovery')

draw_kaplan_meier(scorename=score_name,
                        df=y_test,
                        save_plot=True,
                        add_risk_counts=False,
                        trialname='Validation')


# ## Forest Plots

# In[7]:


draw_forest_plot(time='os.time',
                    event='os.evnt',
                    df=y_train[~y_train['Clinical Trial'].isin(['AML02','AML08','AML05'])],
                    trialname='Discovery',
                    scorename=score_name,
                    save_plot=True)
draw_forest_plot(time='efs.time',
                    event='efs.evnt',
                    df=y_train[~y_train['Clinical Trial'].isin(['AML02','AML08','AML05'])],
                    trialname='Discovery',
                    scorename=score_name,
                    save_plot=True)


# In[8]:


draw_forest_plot(time='os.time',
                    event='os.evnt',
                    df=y_test,
                    trialname='Validation',
                    scorename=score_name,
                    save_plot=True)


# In[9]:


draw_forest_plot(time='efs.time',
                    event='efs.evnt',
                    df=y_test[y_test['Clinical Trial'].isin(['AML02'])],
                    trialname='Validation',
                    scorename=score_name,
                    save_plot=True)


# ## Patient Characteristics Table

# In[10]:


from tableone import TableOne


# In[11]:


df_all = pd.concat([y_train,y_test], join='outer',keys=['StJude (Validation)','COG (Discovery)']).reset_index(level=0, names='cohort')


# In[12]:


columns = ['Age (years)','Age group (years)','Sex','Race or ethnic group',
            'Hispanic or Latino ethnic group', 'MRD 1 Status',
            'Vital Status', 'Leucocyte counts (10⁹/L)',
            'Risk Group','FLT3 ITD']


# In[13]:


mytable_all = TableOne(df_all, columns,overall=False,missing=False,pval=True,
                        order={'FLT3 ITD':['Yes','No'],
                                'Race or ethnic group':['White','Black or African American','Asian'],
                                'MRD 1 Status': ['Positive'],
                                'Risk Group': ['High Risk', 'Standard Risk'],
                                'FLT3 ITD': ['Yes'],
                                'Leucocyte counts (10⁹/L)': ['≥30'],
                                'Age group (years)': ['≥10']},
                        groupby='cohort')
mytable_all.tabulate(tablefmt="html", headers=['',"",'Validation','Discovery', 'p-value'])


# In[14]:


columns2 = ['Age (years)','Age group (years)','Sex','Race or ethnic group',
            'Hispanic or Latino ethnic group', 'MRD 1 Status',
            'Vital Status', 'Leucocyte counts (10⁹/L)',
            'Risk Group','FLT3 ITD','Clinical Trial']

mytable_all = TableOne(y_train, columns2,overall=True,missing=True,pval=False,
                        order={'FLT3 ITD':['Yes','No'],
                                'Race or ethnic group':['White','Black or African American','Asian'],
                                'MRD 1 Status': ['Positive'],
                                'Risk Group': ['High Risk', 'Standard Risk'],
                                'FLT3 ITD': ['Yes'],
                                'Leucocyte counts (10⁹/L)': ['≥30'],
                                'Age group (years)': ['≥10'],
                                'Clinical Trial': ['AAML1031','AAML0531','AML05','AAML03P1']},
                        groupby=None)
mytable_all.tabulate(tablefmt="html", headers=['',"",'Missing','Discovery Cohort'])


# In[15]:


mytable_all = TableOne(y_test, columns2,overall=True,missing=True,pval=False,
                        order={'FLT3 ITD':['Yes','No'],
                                'Race or ethnic group':['White','Black or African American','Asian'],
                                'MRD 1 Status': ['Positive'],
                                'Risk Group': ['High Risk', 'Standard Risk'],
                                'FLT3 ITD': ['Yes'],
                                'Leucocyte counts (10⁹/L)': ['≥30'],
                                'Age group (years)': ['≥10'],
                                'Clinical Trial': ['AML02']},
                        groupby=None)
mytable_all.tabulate(tablefmt="html", headers=['',"",'Missing','Validation Cohort'])


# In[16]:


mytable_all = TableOne(y_train, columns,overall=True,missing=True,pval=False,
                        order={'FLT3 ITD':['Yes','No'],
                                'Race or ethnic group':['White','Black or African American','Asian'],
                                'MRD 1 Status': ['Positive'],
                                'Risk Group': ['High Risk', 'Standard Risk'],
                                'FLT3 ITD': ['Yes'],
                                'Leucocyte counts (10⁹/L)': ['≥30'],
                                'Age group (years)': ['≥10']},
                        groupby='Clinical Trial')
mytable_all.tabulate(tablefmt="html")


# In[ ]:




