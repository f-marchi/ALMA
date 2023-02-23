#!/usr/bin/env python
# coding: utf-8

# # Clinical Data Processing

# ## Where the data at?

# In[1]:


input_path = '../Data/Processed_Data/2_MethylData_Processing_Output.pkl'
clinicaldata_path = '../Data/Raw_Data/Clinical_Data/'
output_path = '../Data/Processed_Data/'


# ## Load Methyl Data

# In[2]:


import pandas as pd

df_methyl = pd.read_pickle(input_path).T.reset_index(level=0, names='Batch')


# ## Add Labels/Clinical Outcome Data

# In[3]:


from FM_Functions.Clinical_Data_CleanUp import *


# In[4]:


labels_cog,labels_aml02,labels_aml08,labels_aml05 = combine_and_index_clinicaldata()


# In[5]:


labels_aml02 = clean_aml02(labels_aml02)
labels_aml08 = clean_aml08(labels_aml08)
labels_cog = clean_cog(labels_cog)
labels_aml05 = clean_aml05(labels_aml05)


# In[6]:


# Combine all clinical data labels
df = pd.concat([labels_aml02,labels_aml08,labels_cog,labels_aml05],axis=0,join='outer')


# In[7]:


def match_methyldata(df):
    """Remove samples that are not in the methyl data"""
    df = df[df.index.isin(df_methyl.index)]
    return df

df = match_methyldata(df)


# In[8]:


# Select only samples from the AML0531 clinical trial (GSE124413)
a = df_methyl[df_methyl['Batch'].isin(['GSE124413_AAML0531'])]
b = df[df.index.isin(a.index)]

# Select control samples from GSE124413
control_0531 = a[~a.index.isin(b.index)]

control_0531['Sample Type'] = 'Bone Marrow Normal'

# Combine control samples with clinical trial samples
df_ = pd.concat([df, control_0531['Sample Type'].to_frame()],axis=0,join='outer')


# ## Remove Samples based on Certain Clinical Features

# ### Remove Relapse Samples

# In[9]:


df1 = df_[~df_['Sample Type'].isin(['Relapse','Recurrent Blood Derived Cancer - Bone Marrow',
                                 'Recurrent Blood Derived Cancer - Peripheral Blood'
                                  ])]
print(f'COG: Of {df_.shape[0]} samples, {df_.shape[0]-df1.shape[0]} matched, yielding {df1.shape[0]} samples after filtering')


# ### Remove Normal Samples

# In[10]:


df2 = df1[~df1['Sample Type'].isin(['Bone Marrow Normal','Blood Derived Normal'])]
print(f'COG: Of {df1.shape[0]} samples, {df1.shape[0]-df2.shape[0]} matched, yielding {df2.shape[0]} samples after filtering')


# ### Remove Duplicate Samples

# In[11]:


df3 = df2[~df2['Patient_ID'].duplicated(keep='last')]
print(f'COG: Of {df2.shape[0]} samples, {df2.shape[0]-df3.shape[0]} matched, yielding {df3.shape[0]} samples after filtering')


# ### Remove Samples from AAML03P1 and CCG2961

# In[12]:


# df4 = df3[df3['Clinical Trial'].isin(['AAML1031','AAML0531','AML02','AML08'])]
# print(f'COG: Of {df3.shape[0]} samples, {df3.shape[0]-df4.shape[0]} matched, yielding {df4.shape[0]} samples after filtering')


# ### Select Control samples

# In[13]:


# controls = df_[df_['Sample Type'].isin(['Bone Marrow Normal'])]

# # Combine control samples with clinical trial samples
# df4 = pd.concat([df3,controls],axis=0,join='outer')

# print(f'COG: Of {df3.shape[0]} samples, {df3.shape[0]-df4.shape[0]} matched, yielding {df4.shape[0]} samples after filtering')


# ## Save Files

# In[14]:


df3['Sample Type'].value_counts(dropna=False)


# In[15]:


df3['Clinical Trial'].value_counts(dropna=False)


# In[16]:


output = df3.join(df_methyl,how='left')


# In[17]:


x = output.iloc[:,df3.shape[1]+1:]
y = output.iloc[:,0:df3.shape[1]+1]


# In[18]:


x.to_pickle(output_path+'x.pkl')
y.to_csv(output_path+'y.csv')


# ## The End
