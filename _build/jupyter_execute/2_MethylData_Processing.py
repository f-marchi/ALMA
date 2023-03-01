#!/usr/bin/env python
# coding: utf-8

# # Methyl Array Data Processing

# ## Data Sources

# Here is a summary of the datasets used:
# 
# 1. AAML1031 ([GSE190931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190931)), __n=1048__
# 2. AAML0531 ([GSE124413](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124413)), __n=500__
# 3. Japanese AML-05 ([GSE133986](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133986)), __n=64__
# 4. AML_TARGET450K from [Genomic Data Commons](https://portal.gdc.cancer.gov/), __n=317__
# 5. AML02,08, liver and control samples from Dr. Jatinder Lamba's Repository, __n=324__

# ## Load Data

# In[1]:


import methylcheck
from pathlib import Path

# AAML0531 Data

filepath1 = Path('../Data/Raw_Data/COG_AAML0531_AAML03P1_GSE124413') # Path to AAML0531 data
betas1 = methylcheck.load(filepath1) # Load AAML0531 data

# AAML1031 Data

filepath2 = Path('../Data/Raw_Data/COG_AAML1031_GSE190931/') # Path to AAML1031 data
betas2 = methylcheck.load(filepath2) # Load AAML1031 data

# Japanese AML05 Data

filepath3 = Path('../Data/Raw_Data/AML05_JapaneseTrial_GSE133986/') # Path to Japanese AML05 data
betas3 = methylcheck.load(filepath3) # Load AML05 data

# TARGET450k data from Genomics Data Commons

filepath4 = Path('../Data/Raw_Data/GDC_TARGET_AML_Methyl450K')
betas4 = methylcheck.load(filepath4) # Load TARGET450k data

# AML02 and 08 450K Data

filepath5 = Path('../Data/Raw_Data/LambaPrivate_StJude_AML02_AML08_Methyl450k') # Path to data
betas5 = methylcheck.load(filepath5) # Load data and metadata


# ## Joining Dataframes

# In[2]:


import pandas as pd

betas = pd.concat([betas1,betas2,betas3,betas4, betas5],
                    keys=['GSE124413_AAML0531','GSE190931_AAML1031',
                    'GSE133986_AML05','GDC_TARGET_AML','StJude_AML02_AML08'],
                    join='inner', axis=1)

print(f' Dataset (df) contains {betas.shape[0]} rows (CpG probes) and {betas.shape[1]} columns (samples).')


# ## Step 1. Filtering Sub-Optimally Designed Probes

# There are several critera for exclusion of probes: Areas that have polymorphisms, cross-hybridization, repeat sequence elements, or base color changes can affect probe quality. 
# 
# Below are publications that have benchmarked probe quality and have provided lists of probes to exclude:
# 
# - [Chen2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3592906/), [Price2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3740789/), [Naeem2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3943510/), [DacaRoszak2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4659175/), [Zhou2016](https://academic.oup.com/nar/article/45/4/e22/2290930)

# We will remove probes based on _Zhou et. al., 2016_. See figure 5 of their paper for detailed description.
# 
# Please see function below `exclude_suboptimal_probes()` for paper information and where to find the annotations.

# In[3]:


def exclude_suboptimal_probes(betas):
    '''This function removes proves listed as sub-optimal according to:
    
    Zhou, W., Laird, P. W. & Shen, H.. Comprehensive characterization,
    annotation and innovative use of Infinium DNA methylation BeadChip probes.
    Nucleic Acids Research gkw967 (2016).
    doi:10.1093/nar/gkw967

    For the .tsv file containing the annotated probes, download the paper's
    supplementary material.
    '''
    zhou2016_probes = pd.read_csv('../Data/UnreliableProbesList_Zhou2016/EPIC.anno.GRCh38.tsv', sep='\t',index_col=0)
    unreliable_probes = list(zhou2016_probes[zhou2016_probes['MASK.general'] == True].index)
    betas_ = betas[~betas.index.isin(unreliable_probes)]
    print(f'Of {betas.shape[0]} probes, {betas.shape[0]-betas_.shape[0]} matched, yielding {betas_.shape[0]} after filtering')
    return(betas_)


# With the function that we created above, we can now remove suboptimal probes:

# In[4]:


df1 = exclude_suboptimal_probes(betas)


# ## Step 2. Filtering sex-linked probes and control probes

# In[5]:


df2 = methylcheck.exclude_sex_control_probes(df1, '450k', no_sex=True, no_control=True, verbose=True)


# ## Step 3. Exclude samples that Illumina QC categorizes as FAIL(pval)

# Here we will:
# 
# 1. Run Illumina's quality control and save reports as excel files
# 2. Load QC reports
# 3. Exclude samples that Illumina categorizes as FAIL(pval) for not meeting the condition: (p ≤ 0.05) > 80% probes

# ### Run Illumina's Quality Control (QC)

# 
# ```methylcheck``` is geared toward quality control of processed data. To this end, there is a helpful function that summarizes performance of control probes ([details on control probes here](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium_hd_methylation/beadarray-controls-reporter-user-guide-1000000004009-00.pdf)). To run this function, the control_probes.pkl file output from ```methylprep``` is required. This report ensures that the chemistry (bisulfite conversion, target specificity, hybridization, staining, etc.) and machine readings are acceptable.
# 
# QC reports are saved as excel sheets in each _filepath_ directory.

# In[6]:


# File path for 1031 has to be uploaded separately since it was processed in batches due to size

filepath2_1 = Path('../Data/Raw_Data/COG_AAML1031_GSE190931/GPL21145_1')
filepath2_2 = Path('../Data/Raw_Data/COG_AAML1031_GSE190931/GPL21145_2')
filepath2_3 = Path('../Data/Raw_Data/COG_AAML1031_GSE190931/GPL21145_3')


# Note: the cell below will generate Illumina's QC report for all datasets and save them in their respective filepaths.
# 
# ```{note}
# It may take >5min to run the cell below. If you don't need the reports in excel format (if they are already there), then skip the cell below.
# ```

# In[7]:


# # AAML0531 QC Report
# methylcheck.controls_report(filepath=filepath1)
# # AAML1031 QC Report
# methylcheck.controls_report(filepath=filepath2_1)
# methylcheck.controls_report(filepath=filepath2_2)
# methylcheck.controls_report(filepath=filepath2_3)
# # AML05 QC Report
# methylcheck.controls_report(filepath=filepath3)
# # TARGET450K QC Report
# methylcheck.controls_report(filepath=filepath4)
# # StJude AML02 and 08 QC Report
# methylcheck.controls_report(filepath=filepath5)


# ### Load QC Reports

# In[8]:


# Load QC reports
qc_table1 = pd.read_excel(str(filepath1)+'/GSE124413_QC_Report.xlsx', index_col=0)
qc_table2_1 = pd.read_excel(str(filepath2_1)+'/GPL21145_1_QC_Report.xlsx', index_col=0)
qc_table2_2 = pd.read_excel(str(filepath2_2)+'/GPL21145_2_QC_Report.xlsx', index_col=0)
qc_table2_3 = pd.read_excel(str(filepath2_3)+'/GPL21145_3_QC_Report.xlsx', index_col=0)
qc_table3 = pd.read_excel(str(filepath3)+'/AML05_Japanese_Trial_QC_Report.xlsx', index_col=0)
qc_table4 = pd.read_excel(str(filepath4)+'/TARGET450k_GDC_QC_Report.xlsx', index_col=0)
qc_table5 = pd.read_excel(str(filepath5)+'/StJude_AML02_AML08_QC_Report.xlsx', index_col=0)

# Merge batches
qc_table = pd.concat([qc_table1.iloc[1:], qc_table2_1.iloc[1:],qc_table2_2.iloc[1:], 
                    qc_table2_3.iloc[1:], qc_table3.iloc[1:], qc_table4.iloc[1:],
                    qc_table5.iloc[1:]], axis=0) # Merge batches


# ### Exclude Samples that Failed Illumina's QC

# Below, you will find the samples that failed QC, along with some description. This column is calculated by checking that all of the QC columns are above a minimum threshold. This threshold is adjustable with the passing argument (set to 0.7 by default). However, this criteria seems to be too unforgiving, so we will disregard it as it seems that core research facilities do not follow it. Another strict metric implemented by Illumina is FAIL by pval, which happens if the pOOBAH (detection p-values) in more than 20% of probes fail. Recall that pOOBAH simply measures how likely it is that a given signal is background fluorescence. There are a few methods of calculating these detection p-values.  ```SeSAMe``` and ```methylprep``` implement this method, known as P value Out Of Band probes for Array Hybridization, where they use the OOB signal of all type I probes to calculate an empirical cumulative distribution function.
# 
# In other words, we will exclude samples that Illumina QC categorizes as FAIL(pval) for __not__ meeting the condition: (pOOBAH ≤ 0.05) > 80% probes

# In[9]:


qc_failed = qc_table[qc_table['Result'].str.contains('pval')][["Result", 'Why Failed']]
qc_failed


# In[10]:


df3 = df2.drop(list(qc_failed.index),level=1, axis=1)
print(f'COG: {df2.shape[1] - df3.shape[1]} sample(s) removed because: (pOOBAH ≤ 0.05) > 80% probes')


# ## Step 4. Exclude CpG probes that contain more than 5% of missing values

# In[11]:


def probe_cutoff(qc_betas, threshold):
    qc_betas2 = qc_betas.dropna(axis=0, thresh = int(threshold*qc_betas.shape[1]))
    print(f'{qc_betas.shape[0] - qc_betas2.shape[0]} probe(s) removed because of >5% missing values')
    return(qc_betas2)


# In[12]:


df4 = probe_cutoff(df3, threshold=0.95)


# ## Step 5. Interpolate remaining missing values linearly

# 
# We still have probes with missing values that are below our 5% threshold. To fix that, we interpolate remaining beta values values linearly.
# - Linear interpolation means that a missing CpG value for a particular sample will be filled with the median of the values of the two adjacent samples. The reason behind this is the high concordance _usually_ seen in the methylation profile of neighboring samples.
# - Inevitably, this step adds arbitration to the data cleaning process.

# In[13]:


df5 = df4.interpolate(axis=0).interpolate(axis=0, limit_direction='backward').round(3)


# Great! In summary:

# In[14]:


print(f' Merged df5 dataset contains {df5.shape[0]}'
+ f' rows (CpG probes) and {df5.shape[1]} columns (samples).')


# ## Step 6. Add Sample Metadata and Clinical Data

# In[15]:


from FM_Functions.Clinical_Data_CleanUp import *
clinical_data = combine_and_index_clinicaldata() # Load and merge all clinical data


# ## Save Files

# In[17]:


df5.to_pickle('../Data/Processed_Data/2_MethylData_Processing_Output.pkl')


# ## End
