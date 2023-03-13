#!/usr/bin/env python
# coding: utf-8

# # Clinical Data Processing

# ## Where the data at?

# In[1]:


input_path = '../Data/Processed_Data/Methyl_Array_Processed/'
clinicaldata_path = '../Data/Raw_Data/Clinical_Data/'
output_path = '../Data/Processed_Data/'


# ## Load Methyl Data

# In[2]:


import pandas as pd

df_methyl = pd.read_pickle(
    input_path+'2_MethylData_Processing_Output.pkl').T.reset_index(level=0, names='Batch')

print(
    f' Dataset (df) contains {df_methyl.shape[1]} rows (mC sites) and {df_methyl.shape[0]} columns (samples).')


# ## Add Labels/Clinical Outcome Data

# In[3]:


# Import functions to clean up clinical data
from FM_Functions.Clinical_Data_CleanUp import *

# Combine all clinical data files into one dataframe and indexes it by the sample ID
labels_cog, labels_aml02, labels_aml08, labels_aml05 = combine_and_index_clinicaldata()

# Clean up and adjust clinical data labels
labels_aml02 = clean_aml02(labels_aml02)
labels_aml08 = clean_aml08(labels_aml08)
labels_cog = clean_cog(labels_cog)
labels_aml05 = clean_aml05(labels_aml05)

# Combine all clinical data labels
df = pd.concat([labels_aml02, labels_aml08, labels_cog,
               labels_aml05], axis=0, join='outer')

# Remove samples that are not in the methyl dataset
df = df.loc[df.index.isin(df_methyl.index)]

# Label control samples from the AML0531 clinical trial (GSE124413) as 'Bone Marrow Normal'


def label_control_samples(df_methyl, df):
    """
    This function labels control samples from the AML0531 clinical trial (GSE124413) as 'Bone Marrow Normal'
    and combines them with the clinical trial samples.
    """
    a = df_methyl[df_methyl['Batch'].isin(['GSE124413_AAML0531'])]
    b = df[df.index.isin(a.index)]
    control_0531 = a[~a.index.isin(b.index)]
    control_0531['Sample Type'] = 'Bone Marrow Normal'
    df_ = pd.concat(
        [df, control_0531['Sample Type'].to_frame()], axis=0, join='outer')
    return df_


df_ = label_control_samples(df_methyl, df)


# ## Remove Samples based on Certain Clinical Features

# ### Remove Relapse Samples

# In[4]:


df1 = df_[~df_['Sample Type'].isin(['Relapse', 'Recurrent Blood Derived Cancer - Bone Marrow',
                                    'Recurrent Blood Derived Cancer - Peripheral Blood'])]

print(
    f'Out of {df_.shape[0]} samples, {df_.shape[0]-df1.shape[0]} matched, yielding {df1.shape[0]} samples after filtering')


# ### Remove Control/Normal Samples

# In[5]:


df2 = df1[~df1['Sample Type'].isin(
    ['Bone Marrow Normal', 'Blood Derived Normal'])]
print(
    f'Out of {df1.shape[0]} samples, {df1.shape[0]-df2.shape[0]} matched, yielding {df2.shape[0]} samples after filtering')


# ### Remove Duplicate Samples

# In[6]:


df3 = df2[~df2['Patient_ID'].duplicated(keep='last')]
print(
    f'Out of {df2.shape[0]} samples, {df2.shape[0]-df3.shape[0]} matched, yielding {df3.shape[0]} samples after filtering')


# ## Save Files

# In[7]:


output = df3.join(df_methyl,how='left') # Join clinical data with methyl data

x = output.iloc[:,df3.shape[1]+1:] # Select only methyl data
y = output.iloc[:,0:df3.shape[1]+1] # Select only clinical data

x.to_pickle(output_path+'x.pkl') # Save methyl data
y.to_csv(output_path+'y.csv') # Save clinical data

print(
    f'Successfuly saved methyl data in x.pkl and clinical data in y.csv.\nPath: {output_path}')


# ### Save Control and Relapse Data Separately

# In[8]:


controls = df_[df_['Sample Type'].isin(['Bone Marrow Normal'])]

relapse = df_[df_['Sample Type'].isin(['Relapse', 'Recurrent Blood Derived Cancer - Bone Marrow',
                                       'Recurrent Blood Derived Cancer - Peripheral Blood'])]

# Merge control and relapse samples
t = pd.concat([controls, relapse], axis=0, join='outer',
              names=['Control', 'Relapse'])

# Join clinical data with methyl data
t2 = df_methyl.join(t, how='right')

# Save merged control and relapse samples
t2.to_pickle(output_path+'control_relapse.pkl')

print(
    f'Successfuly saved {controls.shape[0]} control samples and {relapse.shape[0]} relapse samples.\nPath: {output_path}')


# ## Watermark

# In[9]:


get_ipython().run_line_magic('load_ext', 'watermark')


# In[10]:


# produce a list of the loaded modules
get_ipython().run_line_magic('watermark', '-v -p pandas')

