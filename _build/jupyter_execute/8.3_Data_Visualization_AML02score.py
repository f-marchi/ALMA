#!/usr/bin/env python
# coding: utf-8

# # Data Visualization

# ## Where the data at?

# In[1]:


## Where the data at?
input_path = '../Data/Processed_Data/Cell_Deconvolution/'
output_path = '../Data/Processed_Data/Cell_Deconvolution/'


# ## Which score is this?

# In[2]:


score_name = 'ARIC_mC_score'


# ## Import Libraries

# In[3]:


import pandas as pd
import seaborn as sns

# Set theme
sns.set_theme(style='white')

# Import Plotting Functions
from FM_Functions.Data_Visualization import *


# ## Load Datasets

# In[4]:


import pandas as pd

y_train = pd.read_csv(input_path+'y_plus_cibersortx_ARICresults_'+ score_name +'.csv', index_col=0)
y_test = pd.read_csv(input_path+'y_plus_cibersortx_ARICresults_'+ score_name +'_test.csv', index_col=0)


# In[17]:


y_train['NK Cell Presence'] = pd.cut(x=y_train['NK'],
                                            bins=[-np.inf, 0.001, np.inf],
                                            labels=['NK Cells -', 'NK Cells +'])

y_test['NK Cell Presence'] = pd.cut(x=y_test['NK'],
                                    bins=[-np.inf, 0.001, np.inf],
                                    labels=['NK Cells -', 'NK Cells +'])


# In[21]:


y_test['NK Cells Activated'].describe()


# In[16]:


y_train['NK'].describe()


# In[13]:


y_train['NK Cell Presence'].value_counts(dropna=False)


# ## Kaplan Meiers

# In[22]:


draw_kaplan_meier(scorename=score_name,
                        df=y_train,
                        save_plot=True,
                        add_risk_counts=False,
                        trialname='Discovery')

draw_kaplan_meier(scorename=score_name,
                        df=y_test,
                        save_plot=True,
                        add_risk_counts=False,
                        trialname='Validation')


# ## Forest Plots

# In[6]:


draw_forest_plot(time='os.time',
                    event='os.evnt',
                    df=y_train,
                    trialname='Discovery',
                    scorename=score_name,
                    save_plot=False)
draw_forest_plot(time='efs.time',
                    event='efs.evnt',
                    df=y_train,
                    trialname='Discovery',
                    scorename=score_name,
                    save_plot=False)



# In[7]:


draw_forest_plot(time='os.time',
                    event='os.evnt',
                    df=y_test,
                    trialname='Validation',
                    scorename=score_name,
                    save_plot=True)
draw_forest_plot(time='efs.time',
                    event='efs.evnt',
                    df=y_test,
                    trialname='Validation',
                    scorename=score_name,
                    save_plot=True)


# ## Patient Characteristics Table

# In[8]:


from tableone import TableOne


# In[9]:


df_all = pd.concat([y_train,y_test], join='outer',keys=['StJude (Discovery)','COG (Validation)']).reset_index(level=0, names='cohort')


# In[10]:


columns = ['Age (years)','Age group (years)','Sex','Race or ethnic group',
            'Hispanic or Latino ethnic group', 'MRD 1 Status',
            'Vital Status', 'Leucocyte counts (10⁹/L)',
            'Risk Group','FLT3 ITD']


# In[11]:


mytable_all = TableOne(df_all, columns,overall=False,missing=False,pval=False,
                        order={'FLT3 ITD':['Yes','No'],
                                'Race or ethnic group':['White','Black or African American','Asian'],
                                'MRD 1 Status': ['Positive'],
                                'Risk Group': ['High Risk', 'Standard Risk'],
                                'FLT3 ITD': ['Yes'],
                                'Leucocyte counts (10⁹/L)': ['≥30'],
                                'Age group (years)': ['≥10']},
                        groupby='cohort')
mytable_all.tabulate(tablefmt="html", headers=[score_name,"",'COG (Validation)','StJude (Discovery)'])


# In[12]:


columns2 = ['Age (years)','Age group (years)','Sex','Race or ethnic group',
            'Hispanic or Latino ethnic group', 'MRD 1 Status',
            'Leucocyte counts (10⁹/L)', 'BM Leukemic blasts (%)',
            'Risk Group', 'Clinical Trial','FLT3 ITD']


# In[13]:


mytable_cog = TableOne(y_test, columns2,
                        overall=False, missing=True,
                        pval=True, pval_adjust=False,
                        htest_name=True,dip_test=True,
                        tukey_test=True, normal_test=True,

                        order={'FLT3 ITD':['Yes','No'],
                                'Race or ethnic group':['White','Black or African American','Asian'],
                                'MRD 1 Status': ['Positive'],
                                'Risk Group': ['High Risk', 'Standard Risk'],
                                'FLT3 ITD': ['Yes'],
                                'Leucocyte counts (10⁹/L)': ['≥30'],
                                'Age group (years)': ['≥10']},
                        groupby=score_name + ' Categorical')
mytable_cog.tabulate(tablefmt="html", 
                        headers=[score_name,"",'Missing','High','Low','p-value','Statistical Test'])


# In[14]:


columns3 = ['Age (years)','Age group (years)','Sex','Race or ethnic group',
            'Hispanic or Latino ethnic group', 'MRD 1 Status',
            'Leucocyte counts (10⁹/L)', 'BM Leukemic blasts (%)',
            'Risk Group','FLT3 ITD',
            'Treatment Arm']


# In[15]:


mytable_aml02 = TableOne(y_train, columns3,
                        overall=False, missing=True,
                        pval=True, pval_adjust=False,
                        htest_name=True,dip_test=True,
                        tukey_test=True, normal_test=True,

                        order={'FLT3 ITD':['Yes','No'],
                                'Race or ethnic group':['White','Black or African American','Asian'],
                                'MRD 1 Status': ['Positive'],
                                'Risk Group': ['High Risk', 'Standard Risk'],
                                'FLT3 ITD': ['Yes'],
                                'Leucocyte counts (10⁹/L)': ['≥30'],
                                'Age group (years)': ['≥10']},
                        groupby= score_name + ' Categorical')
mytable_aml02.tabulate(tablefmt="html",
                        headers=[score_name,"",'Missing','High','Low','p-value','Statistical Test'])


# ## Bar Plots

# In[16]:


import numpy as np


# In[21]:


# Set up the matplotlib figure
sns.set_theme(style='white')
f, axs = plt.subplots(2, 1, sharex=True, figsize=(8,7))

# Define plots
sns.histplot(data=y_train,x=score_name, hue=score_name + ' Categorical', ax=axs[0], bins=75)
sns.histplot(data=y_test,x=score_name, hue=score_name + ' Categorical', ax=axs[1], bins=50)

# Set specs
cutoff = np.quantile(y_train[score_name],0.75)

for i in range(2):
    axs[i].axvline(cutoff, linestyle="dotted",color='red', label='Discovery 3ʳᵈ Quartile ('+ round(cutoff,3).astype(str)+ ')')

axs[0].set_title(' Discovery', loc='center', pad=5, fontsize=11)
axs[1].set_title(' Validation', loc='center', pad=5, fontsize=11)

axs[1].legend()
    # Define Plot Specs
plt.subplots_adjust(wspace=0, hspace=0.1)
plt.suptitle(score_name + ' Classifier Strategy: Discovery 3ʳᵈ Quartile ('+ round(cutoff,3).astype(str)+ ')',
                 fontsize='medium', y=0.95,
                 fontweight='bold')
plt.savefig('../Figures/Bar_Plots/'+score_name+'_Classifier_Strategy.png',
                    bbox_inches='tight', dpi=300)

plt.show()


# ## Box Plots

# In[24]:


draw_boxplot(df=y_train,x='Risk Group', y=score_name,
                order=['High Risk', 'Standard Risk', 'Low Risk'],
                trialname='Discovery', hue=score_name +' Categorical',
                save_plot=True, figsize=None)

draw_boxplot(df=y_test,x='Risk Group', y=score_name,
                order=['High Risk', 'Standard Risk', 'Low Risk'],
                trialname='Validation', hue=score_name + ' Categorical',
                save_plot=True, figsize=None)


# In[25]:


draw_boxplot(df=y_train,x='MRD 1 Status', y=score_name,
                order=['Positive','Negative'],
                trialname='Discovery', hue=score_name + ' Categorical',
                save_plot=True, figsize=None)

draw_boxplot(df=y_test,x='MRD 1 Status', y=score_name,
                order=['Positive','Negative'],
                trialname='Validation', hue=score_name + ' Categorical',
                save_plot=True, figsize=None)


# In[26]:


draw_boxplot(df=y_train,x='Primary Cytogenetic Code', y=score_name,
                order='auto',
                trialname='Discovery', hue=score_name + ' Categorical',
                save_plot=True, figsize=None)

draw_boxplot(df=y_test,x='Primary Cytogenetic Code', y=score_name,
                order='auto',
                trialname='Validation', hue=score_name + ' Categorical',
                save_plot=True, figsize=None)


# ## End
