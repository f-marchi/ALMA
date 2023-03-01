#!/usr/bin/env python
# coding: utf-8

# # Immune Cell Deconvolution

# ## Where the data at?

# In[3]:


input_path = '../Data/Processed_Data/'
output_path = '../Data/Processed_Data/Blood_Deconvolution_ARIC/'


# ## Load AML Dataset

# In[4]:


import pandas as pd

x = pd.read_pickle(input_path+'x.pkl')
y = pd.read_csv(input_path+'y.csv', index_col=0)


# ## Train-Test Split

# To avoid data leakage and maximize the independence of the validation cohort (test dataset), we will split the data by clinical trial. The validation cohort will be St. Jude Children's led trials (AML02 and AML08), and the training cohort will be all other trials.

# In[5]:


# Split data into training and test sets by clinical trial
y_train = y[~y['Clinical Trial'].isin(['AML02','AML08'])]
y_test = y[y['Clinical Trial'].isin(['AML02','AML08'])]

# Select samples in x that are in y_train
x_train = x.loc[y_train.index]
x_test = x.loc[y_test.index]

y['Clinical Trial'].value_counts(dropna=False)


# ## Batch Correction with pyCombat

# - __pyCombat__: a Python tool for batch effects correction in high-throughput molecular data using empirical Bayes methods
# 
# - __Website__: [https://epigenelabs.github.io/pyComBat/](https://epigenelabs.github.io/pyComBat/)
# 
# - __Paper__: [bioRxiv](https://doi.org/10.1101/2020.03.17.995431)

# In[6]:


from combat.pycombat import pycombat
data_corrected = pycombat(x_train.T,y_train['Batch'])
x_train2 = data_corrected.T


# ## Load Reference Dataset

# - __FlowSorted.Blood.EPIC__: An optimized library for reference-based deconvolution of whole-blood biospecimens, __n=49__
# 
# - __GEO__: [GSE110554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554)
# 
# - __PMID__: [29843789](https://www.ncbi.nlm.nih.gov/pubmed/29843789)
# 
# - __Description__:  Bisulphite converted DNA from neutrophils (Neu, n=6), monocytes (Mono, n=6), B-lymphocytes (Bcells, n=6), CD4+ T-cells (CD4T, n=7, six samples and one technical replicate), CD8+ T-cells (CD8T, n=6), Natural Killer cells (NK, n=6), and 12 DNA artificial mixtures (labeled as MIX in the dataset) were hybridised to the Illumina Infinium HumanMethylationEPIC Beadchip v1.0_B4
# 
# - __CSV file__: [Download](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-018-1448-7/MediaObjects/13059_2018_1448_MOESM4_ESM.csv)

# In[7]:


ref = pd.read_csv('https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-018-1448-7/MediaObjects/13059_2018_1448_MOESM4_ESM.csv',
                  index_col=0, skiprows=1)[['CD8T','CD4T','NK','Bcell','Mono','Neu']]

# File has also been downloaded locally under ".../Data/Blood_Reference_PMID29843789" as backup in case Springer link is down

mix = x_train.T
merge = ref.join(mix, how='inner')

# update ref and mix with merge index
ref = ref.loc[merge.index]
mix = mix.loc[merge.index]

# save ref and mix to csv
ref.to_csv(output_path+'ref.csv')
mix.to_csv(output_path+'mix.csv')


# ## Immune Cell Deconvolution with ARIC

# - __ARIC__: Accurate and robust inference of cell type proportions from bulk gene expression or DNA methylation data
# 
# - __Website__: [xwanglabthu.github.io/ARIC/](xwanglabthu.github.io/ARIC/)
# 
# - __PMID__: [34472588](https://pubmed.ncbi.nlm.nih.gov/34472588/)
# 
# - __External Validation__: [A systematic assessment of cell type deconvolution algorithms for DNA methylation data](https://doi.org/10.1093/bib/bbac449)

# In[8]:


from ARIC import *

ARIC(mix_path=output_path+'mix.csv', ref_path=output_path+'ref.csv',
     is_methylation=True, unknown=False)


# In[16]:


# Read deconvolution results

deconv = pd.read_csv(output_path+'mix_prop.csv', index_col=0)

# Get basic statistics of deconvolution results
round(deconv.T.describe(),2)

