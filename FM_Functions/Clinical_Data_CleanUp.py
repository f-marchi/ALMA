"""
This module cleans and adjusts clinical data files from St. Jude and COG AML trials.

"""

import pandas as pd
import numpy as np

__author__ = 'Francisco Marchi, Lamba Lab, University of Florida'
__email__ = 'flourenco@ufl.edu'


clinical_data_path='../Data/Raw_Data/Clinical_Data/'

def merge_index_1031(filepath1 = '/TARGET/TARGET-AML/TARGET_AML_ClinicalData_AML1031_20221108.xlsx',
                     filepath2 = '../Data/Raw_Data/Methyl_Array_EPIC/GSE190931/sample_sheet_meta_data.pkl'):

    # Load clinical data files
    labels_1031 = pd.read_excel(clinical_data_path + filepath1)
    meta        = pd.read_pickle(filepath2)

    # Extract last term of `TARGET USI` by splitting on `-`
    labels_1031['Patient_ID'] = labels_1031['TARGET USI'].str.split('-').str[2]

    # extract patient name from `description` column by splitting on `\n` and then on `_`
    meta['Patient_ID'] = meta['description'].str.split('\n').str[-1].str.split('_').str[-2]

    # Set index to `Patient_ID` and selected only relevant columns
    meta = meta.set_index('Patient_ID').iloc[:,:-1][['fusion','timepoint','Sample_ID']]

    # Set index to `Patient_ID`
    labels_1031 = labels_1031.set_index('Patient_ID')

    # Join the two dataframes
    labels_1031 = labels_1031.join(meta, how='right').set_index('Sample_ID')

    return labels_1031

def merge_index_0531(dir = '../Data/Raw_Data/Clinical_Data/TARGET/TARGET-AML/',
                    filepath1 = 'TARGET_AML_ClinicalData_Discovery_20221108.xlsx',
                    filepath2 = 'TARGET_AML_ClinicalData_Validation_20221108.xlsx',
                    filepath3 = 'TARGET_AML_ClinicalData_LowDepthRNAseq_20221108.xlsx',
                    filepath4 = '../Data/Raw_Data/Methyl_Array_EPIC/GSE124413/sample_sheet_meta_data.pkl',
                    filepath5 = '../Data/Raw_Data/Methyl_Array_EPIC/GSE124413/GSE124413_series_matrix.csv'):
    
    # Load all clinical data files for 0531
    labels_0531_1 = pd.read_excel(dir + filepath1, index_col=0)
    labels_0531_2 = pd.read_excel(dir + filepath2, index_col=0)
    labels_0531_3 = pd.read_excel(dir + filepath3, index_col=0)
    meta          = pd.read_pickle(filepath4)
    meta_matrix   = pd.read_csv(filepath5)

    # Concatenate the two dataframes
    labels_0531 = pd.concat([labels_0531_1, labels_0531_2], axis=0, join='outer').reset_index()

    def remove_the_duplicate_samples_with_more_nulls(df=labels_0531):
        ''' 
        This function removes duplicate samples from the dataframe, keeping the row with fewer NaNs (null values).
        '''    
        # Adding a new column 'nan_count' which is the count of NaNs in each row
        df = df.replace({'NA': np.nan, 'unknown': np.nan, 'Unknown': np.nan})
        df['nan_count'] = df.isnull().sum(axis=1)

        # Sort by 'nan_count' so that rows with fewer NaNs come first
        df = df.sort_values('nan_count')

        # Drop duplicates, keeping the first one (with fewer NaNs)
        df = df.drop_duplicates(subset= 'TARGET USI', keep='first')

        # Remove the 'nan_count' column
        df = df.drop(columns='nan_count')
        return df

    def clean_meta(meta_matrix, meta):
        # Transpose the dataframe and reset index
        transposed_meta = meta_matrix.T.reset_index()

        # Split the index and join to dataframe
        transposed_meta['new_index'] = transposed_meta['index'].str.rsplit(" ", n=1, expand=True)[1]
        
        # Set new header
        transposed_meta.columns = transposed_meta.iloc[0]
        transposed_meta = transposed_meta.drop(transposed_meta.index[0])

        # Rename columns and set index
        transposed_meta = transposed_meta.rename(columns={'!Sample_geo_accession':'GSM_ID', None: 'Patient_ID'}).set_index('GSM_ID')

        # Join with meta DataFrame, select columns and reset index
        meta_cleaned = meta.set_index('GSM_ID').join(transposed_meta)[['Patient_ID', '!Sample_characteristics_ch1','Sample_ID']].reset_index().set_index('Patient_ID')

        # Rename columns
        meta_cleaned.columns = ['GSM_ID', 'Sample Type', 'age','sex', 'Sample_ID']

        return meta_cleaned

    # Adjust and clean metadata
    meta_cleaned = clean_meta(meta_matrix, meta)

    # Remove duplicate samples by keeping the row with fewer NaNs
    labels_0531 = remove_the_duplicate_samples_with_more_nulls()

    # Set index to `TARGET USI` and join with `Gene Fusion.1` column
    labels_0531 = labels_0531.set_index('TARGET USI').join(labels_0531_3[['Gene Fusion.1']], how='outer')

    # extract the last two words of `TARGET USI` by splitting on `-` 
    labels_0531['Tumor Code'] = labels_0531.index.str.split('-').str[1]
    labels_0531['Patient_ID'] = labels_0531.index.str.split('-').str[2]

    # drop duplicates in `Patient_ID` column
    labels_0531 = labels_0531.drop_duplicates(subset='Patient_ID', keep='first').set_index('Patient_ID')

    # join with meta_cleaned dataframe
    labels_0531 = labels_0531.join(meta_cleaned, how='right')

    # drop columns that are not needed
    labels_0531 = labels_0531.drop(columns=['age', 'sex', 'GSM_ID'])

    # Rename values in `Sample Type` column
    labels_0531['Sample Type'] = labels_0531['Sample Type'].replace({'group: tumor': 'Diagnosis',
                                                                    'group: normal': 'Bone Marrow Normal'})

    # Set index to `Sample_ID` to match methylation samples
    labels_0531 = labels_0531.set_index('Sample_ID')

    return labels_0531

# COG/TARGET-AML
def merge_index_cog():
    labels_0531 = pd.read_csv(clinical_data_path + 'MarchiF_ClinicalData_AML0531_03P1.csv',
                                    index_col=67)
    labels_0531['Sample Type'] = 'Diagnosis'
    labels2_1 = pd.read_csv(clinical_data_path + 'MarchiF_ClinicalData_AAML1031.csv',
                                index_col=67).rename(columns={'Timepoint': 'Sample Type'})
    labels2_2 = pd.read_excel(clinical_data_path + 'COGAAML1031.v3.11.9.21.xlsx',
                                    index_col=0)
    labels_1031 = labels2_1.join(labels2_2[['Treatment Arm', 'KRAS', 'Protocol risk group classification']],
                                how='left', on='Patient_ID')
    labels3_1 = pd.read_csv(
                                clinical_data_path + 'clinical_info_gdc_2022-06-28.tsv', sep='\t')
    COG_clinicaldata_all = pd.read_excel(
                                clinical_data_path + 'COG_raw_clinicaldata.xlsx', index_col=1)

    def clean_TARGET450k(labels3, COG_clinicaldata_all):
        l1 = labels3['File Name'].str.rsplit('_', n=1, expand=True).rename({
                0: 'IlmnID'}, axis='columns')
        l2 = l1.join(labels3)
        l3 = l2[l2[1] == 'Grn.idat'].rename(
                {'Case ID': 'TARGET USI'}, axis='columns')
        l4 = l3[~l3['Sample Type'].isin(
                ['Fibroblasts from Bone Marrow Normal'])]
        l5 = l4.set_index('TARGET USI').join(COG_clinicaldata_all, how='left')
        l6 = l5.drop(columns=[1, 'File ID', 'Data Type', 'Project ID',
                                'Data Category', "Unnamed: 0",
                                'File Name']).drop_duplicates(subset=['IlmnID'], keep='last').set_index('IlmnID')
        return (l6)

    labels_gdc_target = clean_TARGET450k(labels3_1, COG_clinicaldata_all)

    labels_cog = pd.concat(
            [labels_0531, labels_1031, labels_gdc_target], axis=0, join='outer')
        
    return (labels_cog)
    
# AML05
def merge_index_aml05():
    labels_aml05 = pd.read_pickle(
    clinical_data_path + 'AML05_sample_sheet_meta_data.pkl'
        ).iloc[:,:-1].set_index('Sample_ID')
    return (labels_aml05)

# AML02
def merge_index_aml02():

        labels5_1 = pd.read_excel(clinical_data_path + 'JLamba-AML02-data-2019-02-12_JKL.Francisco.xlsx',
                              index_col=94).drop(columns=['DONOTUSE: ACCESSION_NBR'])  # main clinical data file
        labels5_2 = pd.read_excel(clinical_data_path + 'AML02_Clinical_Data_Mastersheet.xlsx',
                              index_col=0)  # only used to retrieve pLSC6 values
        labels5_3 = pd.read_excel(clinical_data_path + 'AML02methylation_Lamba_July25.2014Summary_FM_Cleaned.xlsx',
                              index_col=1)  # used to merged clinical data with methyl
        # Join contents of two columns into one

        labels5_3['Sample'] = labels5_3['Sentrix Barcode'].astype(
            str) + '_' + labels5_3['Sample Section'].astype(str)
        labels5_4 = labels5_3[labels5_3['type'] ==
                              'AML02'].reset_index()
        # Replace - with _ in Sample ID

        labels5_4['Sample ID'] = labels5_4['Sample ID'].str.replace(
            '-', '_')
        # Set index to Sample ID

        labels5_4 = labels5_4.set_index('Sample ID')

        # Combine all clinical data files together

        labels5_5 = labels5_1.join(labels5_2[['pLSC6_gb', 'pLSC6_Score', 'DNMT3B']],
                                   how='left', on='U133A.Dx.ID').join(labels5_4['Sample'],
                                                                      how='inner', on='JL.id.x').reset_index(
        ).set_index('Sample')

        # Remove samples that Dr. Lamba requested to be removed

        labels5_6 = labels5_5[labels5_5['Meth450K.Array.ID'].isin(
            ['AML02_119_P022', 'AML02_147_P053', 'AML02_197_P022',
             'AML02_23_P022', 'AML02_39_P022', 'AML02_46_P022',
             'AML02_13_P022']).astype(int) == 0]
        return (labels5_6)

# AML08
def merge_index_aml08():
        labels6_1 = pd.read_excel(clinical_data_path + 'AML08-clinical-AEs-IDs-2022-06-01.xlsx',
                                sheet_name=[0, 2],
                                index_col=0)

        labels6_2 = pd.read_csv(clinical_data_path + 'AML02.AML08.PC.clin.rand.merge.N400.csv',
                                index_col=1)
        # Merge cleaned clinical data files with methylation data file
        labels6_1[0] = labels6_1[0].join(labels6_2['AGE'], how='left')
        labels_aml08 = labels6_1[2]['Meth.Sample.ID'].to_frame().join(labels6_1[0],
                                                                    how='inner').reset_index(
        ).set_index('Meth.Sample.ID')
        labels_aml08.columns = labels_aml08.columns.str.title()
        return (labels_aml08)

# BeatAML
def merge_index_beataml ():
        meta = pd.read_pickle('../Data/Raw_Data/Methyl_Array_EPIC/GSE159907/sample_sheet_meta_data.pkl').iloc[:,:-1]

        # Create a new column with only the content inside [] from column 'Sample_Name'
        meta['LLS_SampleID'] = meta['Sample_Name'].str.extract(r"\[(.*?)\]", expand=False)

        # Set the index to the new column
        meta1 = meta[['tissue','disease_state','LLS_SampleID','Sample_ID']].set_index('LLS_SampleID')

        # Read in the clinical data
        meta2 = pd.read_excel('../Data/Raw_Data/Clinical_Data/BeatAML/BEAT_AML_Raw clinical data_702.Samples.Vizome.xlsx', index_col=3)

        # Join the two dataframes
        labels_beataml = meta1.join(meta2, how='left').reset_index().set_index('Sample_ID')

        return labels_beataml

# AML-TCGA
def merge_index_amltcga():

        # load clinical data from GDC
        clinical_tsv = pd.read_csv('../Data/Raw_Data/Methyl_Array_450k/GDC_TCGA-AML/clinical.tsv', 
                        sep='\t', index_col=0)[['case_submitter_id']].drop_duplicates()

        # extract last 4 digits from case_id to get TCGA Patient ID
        clinical_tsv['TCGA Patient ID'] = clinical_tsv['case_submitter_id'].str[-4:]

        # set index to TCGA Patient ID
        clinical_tsv = clinical_tsv.reset_index().set_index('TCGA Patient ID').sort_index()

        # load meta data from NEJM 2013 paper
        meta = pd.read_excel('../Data/Raw_Data/Clinical_Data/TCGA_LAML/SuppTable01_NEJM2013_TCGA_AML.Paper_Mutation data.xlsx',
                            index_col=1).iloc[1:,:].sort_index()

        # make meta index integers
        meta.index = meta.index.astype(int)
        clinical_tsv.index = clinical_tsv.index.astype(int)

        # join clinical_tsv and meta
        labels_amltcga = clinical_tsv.join(meta, how='left')

        # add `_noid` to end of case_id to match methylation samples
        labels_amltcga['case_id'] = labels_amltcga['case_id']+'_noid'

        # set index to case_id
        labels_amltcga = labels_amltcga.reset_index().set_index('case_id')
        
        return labels_amltcga

# Nordic ALL
def merge_index_nordic_all(): 
        # Load meta data from GSE49031
        meta = pd.read_pickle('../Data/Raw_Data/Methyl_Array_450k/GSE49031/sample_sheet_meta_data.pkl')\
                                .iloc[:,:-1].set_index('Sample_ID')

        # split meta `title` column by the last word
        meta['title'] = meta['title'].str.split().str[-1]

        # Set index to `title`
        meta = meta.reset_index().set_index('title')

        # Load clinical data from paper
        paper = pd.read_excel('../Data/Raw_Data/Clinical_Data/Nordic_ALL/PMID_25729447_Supp_Clinical_Data.xlsx',
                            index_col=0,header=2, sheet_name='Table S7- Verification summary')[['Karyotyping at diagnosisc']]

        # Join meta and paper
        meta = meta.join(paper)

        # Reset index to `Sample_ID`
        meta = meta.reset_index().set_index('Sample_ID')
        return meta

# MDS_tAML
def merge_index_mds_taml():
        meta = pd.read_pickle('../Data/Raw_Data/Methyl_Array_450k/GSE152710/sample_sheet_meta_data.pkl')\
            .iloc[:,:-1].set_index('Sample_ID')
        return meta
    
# Tcell_ALL_GRAAL
def merge_index_all_graal():
        meta = pd.read_pickle('../Data/Raw_Data/Methyl_Array_EPIC/GSE147667/sample_sheet_meta_data.pkl')\
            .iloc[:,:-1].set_index('Sample_ID')
        return meta
        
# GDC_TARGET-ALL
# The clinical data for this dataset does not match the methylation data.
# It is likely that controlled TARGET access is required to access the correct clinical data.
def merge_index_target_all():

        # # Load clinical data from GDC
        # json_clinical_demographic = pd.read_json('../Data/Raw_Data/Methyl_Array_EPIC/GDC_TARGET-ALL/clinical.cases_selection.2023-05-12.json',
        #                             orient='values')
        # # flatten json
        # json_clinical_demographic = pd.json_normalize(json_clinical_demographic['demographic'].dropna())

        # # extract the second to last term from the `submitter_id` column
        # json_clinical_demographic['submitter_id'] = json_clinical_demographic['submitter_id'].str.split('-').str[-1]

        # # extract the first term from the `submitter_id` column by `_`
        # json_clinical_demographic['submitter_id'] = json_clinical_demographic['submitter_id'].str.split('_').str[0]

        # # change `submitter_id` column name to `Patient_ID`
        # json_clinical_demographic = json_clinical_demographic.rename(columns={'submitter_id':'Patient_ID'})

        # # Set index to `submitter_id`
        # json_clinical_demographic = json_clinical_demographic.set_index('demographic_id')['Patient_ID']

        # # Load clinical data from GDC
        # clinical_tsv = pd.read_csv('../Data/Raw_Data/Methyl_Array_EPIC/GDC_TARGET-ALL/clinical.tsv', 
        #                             sep='\t', index_col=0)

        # # Extract the last word from the `case_submitter_id` column by splitting by `-`
        # clinical_tsv['Patient_ID'] = clinical_tsv['case_submitter_id'].str.split('-').str[-1]

        # clinical_tsv = clinical_tsv['Patient_ID']

        # # concat clinical_tsv and json_clinical_demographic
        # clinical = pd.concat([clinical_tsv, json_clinical_demographic], axis=0, join='outer')

        # # Set index to `Patient_ID`
        # clinical = clinical.reset_index().set_index('Patient_ID')

        # # Load clinical data from paper
        # paper = pd.read_excel('../Data/Raw_Data/Clinical_Data/ALL_P3_TARGET/41586_2018_436_MOESM4_ESM.xlsx',
        #                     sheet_name='ST2 Cohort', index_col=0)

        # # # Join clinical data from paper and GDC
        # labels_alltarget = clinical.join(paper, how='right')

        meta = pd.read_pickle('../Data/Raw_Data/Methyl_Array_EPIC/GDC_TARGET-ALL/sample_sheet_meta_data.pkl')\
                              .set_index('Sample_ID')['Sentrix_ID'].to_frame()
        
        return meta


def clean_aml02(df):
    """
    This module cleans and adjusts clinical data files from St. Jude AML02 trial.

    Parameters:
    -----------
    df:object
        AML02 trial clinical outcome dataframe.

    Returns:
    --------
        Cleaned and adjusted AML02 dataframe.

    """
    df = df.rename(columns={'efs.evst': 'efs.evnt', 'os.evst': 'os.evnt', 'efs.evtm': 'efs.time',
                            'os.evtm': 'os.time', 'MRD22': 'MRD 1 Status', 'initialRisk': 'Risk Group',
                            'CYTres': 'Primary Cytogenetic Code', 'FLT3': 'FLT3 Status', 'STUDY': 'Clinical Trial',
                            'WBC': 'WBC Count (10⁹/L)', 'age.dx': 'Age (years)', 'ARM': 'Study Arm',
                            'SEX': 'Sex', 'RACE': 'Race', 'ETHNICITY_CD': 'Ethnicity',
                            'DXCBCBLASTS': 'Peripheral blasts (%)', 'JL.id.x': 'Patient_ID',
                            'DXBMBLAST': 'BM Leukemic blasts (%)', 'REASONOFFSTUDY_02': 'Reasonoffstudy',
                            'CAUSEOFDEATH': 'Causeofdeath', 'dx02': 'Primarydiagnosis',
                            'FULLKARYOTYPE': 'Karyotype'})  # Note that patient ID column is JL.id.x

    df['FLT3 Status'] = df['FLT3 Status'].replace({'Wild type': 'Wild Type'})
    df['Ethnicity'] = df['Ethnicity'].replace({'Non Hispanic or Latino': 'Not Hispanic or Latino',
                                               'Hispanic or Latino': 'Hispanic or Latino'})
    df['relapse.evnt'] = df['fail.type'].isin(['Relapse']).astype(int)
    df['rel_indfail.evnt'] = df['fail.type'].isin(
        ['Relapse', 'Induction Failure']).astype(int)
    df['Treatment Arm'] = df['Study Arm'].replace(
        {'LDAC': 'Arm A', 'HDAC': 'Arm B', 'Not Randomized': np.nan})
    df['Sex'] = df['Sex'].replace({1: 'Male', 2: 'Female'})
    df['Vital Status'] = df['os.evnt'].replace({0: 'Alive', 1: 'Dead'})
    df = df.replace(
        {'Unknown': np.nan, 'Unkown': np.nan, 'Inevaluable': np.nan})
    df['FLT3 ITD'] = df['FLT3 Status'].replace(
        {'Wild Type': 'No', 'Mutation': 'No', 'ITD': 'Yes'})
    df['Age group (years)'] = pd.cut(x=df['Age (years)'],
                                     bins=[-np.inf, 10, np.inf],
                                     labels=['<10', '≥10'])
    df['Leucocyte counts (10⁹/L)'] = pd.cut(x=df['WBC Count (10⁹/L)'],
                                            bins=[-np.inf, 30, np.inf],
                                            labels=['<30', '≥30'])
    df['Race'] = df['Race'].replace({'Black': 'Black or African American',
                                     'Pacific Islander': 'Native Hawaiian or other Pacific Islander',
                                     'Multiple Race  (NOS)': 'Other'})
    df['Risk Group'] = df['Risk Group'].replace({'High': 'High Risk',
                                                 'Standard': 'Standard Risk',
                                                 'Low': 'Low Risk'})
    df = df.rename(columns={'Race': 'Race or ethnic group',
                            'Ethnicity': 'Hispanic or Latino ethnic group'})
    df['Tissue Type'] = 'Bone Marrow'
    df['Sample Type'] = 'Diagnosis'

    return (df)


def clean_aml08(df):
    """
    This module cleans and adjusts clinical data files from St. Jude AML08 trial.

    Parameters:
    -----------
    df:object
        AML08 trial clinical outcome dataframe.

    Returns:
    --------
        Cleaned and adjusted AML08 dataframe.
    """

    df = df.rename(columns={'Efscensor': 'efs.evnt', 'Oscensor': 'os.evnt',
                            'Ostime': 'os.time.days', 'Efstime': 'efs.time.days',
                            'Ethnicity_Cd': 'Ethnicity', 'Dxbmblast': 'BM Leukemic blasts (%)',
                            'Dxcbcblast': 'Peripheral blasts (%)', 'Gender': 'Sex',
                            'Age': 'Age (years)', 'Wbc': 'WBC Count (10⁹/L)',
                            'Flt3': 'FLT3 Status', 'Aml08.Sp.Id': 'Patient_ID',
                            'Fullkaryotypetesting': 'Karyotype'})
    df['Ethnicity'] = df['Ethnicity'].replace({'Non Spanish speaking, Non Hispanic': 'Not Hispanic or Latino',
                                               'NOS Spanish,Hispanic,Latino': 'Hispanic or Latino'})
    df['os.time'] = df['os.time.days'].astype(float).apply(lambda x: x/365)
    df['efs.time'] = df['efs.time.days'].astype(float).apply(lambda x: x/365)
    df['Clinical Trial'] = 'AML08'
    df['Risk Group'] = df['Dx.Risk'].replace(
        {1: 'Low', 2: 'Standard', 3: 'High'})
    df['MRD 1 Status'] = df['Mrdgrp1'].replace(
        {1.0: 'Negative', 2.0: 'Positive'})
    df['Treatment Arm'] = df['Randomizedarm'].replace({'Clo/AraC': 'Arm A',
                                                       'HD-ADE': 'Arm B',
                                                       'Not Randomized (give HD-ADE)': 'Arm B'})
    df['Vital Status'] = df['os.evnt'].replace({0: 'Alive', 1: 'Dead'})
    df = df.replace(
        {'Unknown': np.nan, 'Unkown': np.nan, 'Inevaluable': np.nan})
    df['FLT3 ITD'] = df['FLT3 Status'].replace(
        {'Wild Type': 'No', 'Mutation': 'No', 'ITD': 'Yes'})
    df['Age group (years)'] = pd.cut(x=df['Age (years)'],
                                     bins=[-np.inf, 10, np.inf],
                                     labels=['<10', '≥10'])
    df['Leucocyte counts (10⁹/L)'] = pd.cut(x=df['WBC Count (10⁹/L)'],
                                            bins=[-np.inf, 30, np.inf],
                                            labels=['<30', '≥30'])
    df['Race'] = df['Race'].replace({'Black': 'Black or African American',
                                     'Pacific Islander': 'Native Hawaiian or other Pacific Islander',
                                     'Multiple Race  (NOS)': 'Other'})
    df['Risk Group'] = df['Risk Group'].replace({'High': 'High Risk',
                                                 'Standard': 'Standard Risk',
                                                 'Low': 'Low Risk'})
    df = df.rename(columns={'Race': 'Race or ethnic group',
                            'Ethnicity': 'Hispanic or Latino ethnic group'})

    df['Tissue Type'] = 'Bone Marrow'
    df['Sample Type'] = 'Diagnosis'

    return (df)


def clean_cog(df):
    """
    This module cleans and adjusts clinical data files from Children's Oncology Group AML trials.

    Parameters:
    -----------
    df:object
        COG/TARGET clinical outcome dataframe.

    Returns:
    --------
        Cleaned and adjusted dataframe.

    """

    df = df.replace({'Unknown': np.nan, 'YES': 'Yes', 'NO': 'No', 'male': 'Male',
                    'female': 'Female', 'bone Marrow': 'Bone Marrow'})
    df['Risk Group'] = df['Risk group'].replace({'10': np.nan, '30': np.nan})
    df['First Event'] = df['First Event'].replace({'Induction failure': 'Induction Failure',
                                                   'Death without remission': 'Death without Remission'})
    df['Treatment Arm'] = df['Gemtuzumab ozogamicin treatment'].replace({'Gentuzumab ozogamicin treatment': 'Arm B',
                                                                         'NO Gentuzumab ozogamicin treatment': 'Arm A'})
    df = df.rename(columns={'FLT3/ITD positive?': 'FLT3 ITD', 'Gender': 'Sex', 'Race': 'Race or ethnic group',
                            'Ethnicity': 'Hispanic or Latino ethnic group',
                            'Bone marrow leukemic blast percentage (%)': 'BM Leukemic blasts (%)',
                            'Protocol': 'Clinical Trial', 'WBC at Diagnosis': 'WBC Count (10⁹/L)',
                            'ISCN': 'Karyotype', 'FAB Category': 'FAB'})
    df['CNS disease'] = df[~df['Clinical Trial'].isin(
        ['AAML1031'])]['CNS disease']
    df['Karyotype Complexity 3'] = pd.to_numeric(df['Cytogenetic Complexity'].replace({'>6': 6, 'More than 6': 6, '4 or more': 4}),
                                                 errors='ignore').replace({range(3, 100): '3 or more'})
    df['Karyotype Complexity 4'] = pd.to_numeric(df['Cytogenetic Complexity'].replace({'>6': 6, 'More than 6': 6, '4 or more': 4}),
                                                 errors='ignore').replace({range(4, 100): '4 or more'})
    df['Karyotype Complexity 6'] = pd.to_numeric(df['Cytogenetic Complexity'].replace({'>6': 6, 'More than 6': 6, '4 or more': 4}),
                                                 errors='ignore').replace({range(6, 100): '6 or more'})
    df['Complex Karyotype'] = df[df['Karyotype Complexity 4'].isin(
        ['4 or more'])]['Karyotype Complexity 4']
    df['Complex Karyotype'] = df['Complex Karyotype'].fillna('Less than 4')
    df['Hispanic or Latino ethnic group'] = df['Hispanic or Latino ethnic group'].replace({
        'Not Hispanic or Latino': 'Not Hispanic or Latino',
        'Hispanic or Latino': 'Hispanic or Latino'})
    df['Age (years)'] = df["Age at Diagnosis in Days"].astype(
        float).apply(lambda x: x/365)
    df['efs.time'] = df['Event Free Survival Time in Days'].astype(
        float).apply(lambda x: x/365)
    df['os.time'] = df['Overall Survival Time in Days'].astype(
        float).apply(lambda x: x/365)
    df['os.evnt'] = df['Vital Status'].map({'Alive': 0, 'Dead': 1})
    df['efs.evnt'] = df['First Event'].isin(
        ['Relapse', 'Induction Failure', 'Death', 'Death without Remission']).astype(int)
    df['treat.arm'] = df['Treatment Arm'].map({'Arm A': 0, 'Arm B': 1})
    df['relapse.evnt'] = df['First Event'].isin(['Relapse']).astype(int)
    df['rel_indfail.evnt'] = df['First Event'].isin(
        ['Relapse', 'Induction Failure']).astype(int)
    df['MRD 1 Status'] = df['MRD at end of course 1'].replace(
        {'Yes': 'Positive', 'No': 'Negative'})
    df['Age group (years)'] = pd.cut(x=df['Age (years)'],
                                     bins=[-np.inf, 10, np.inf],
                                     labels=['<10', '≥10'])
    df['Leucocyte counts (10⁹/L)'] = pd.cut(x=df['WBC Count (10⁹/L)'],
                                            bins=[-np.inf, 30, np.inf],
                                            labels=['<30', '≥30'])
    # df['Tissue Code'] = df['Tissue Code'].replace({'09A': '09A - Primary blood derived cancer - bone marrow',
    #                                                '03A': '03A - Primary blood derived cancer - peripheral blood',
    #                                                '14A': '14A - Bone marrow normal',
    #                                                '04A': '04A - Recurrent blood derived cancer – bone marrow',
    #                                                '41A': '41A - Blood derived cancer- bone marrow, post-treatment',
    #                                                '50A': '50A - Cell line from patient tumor',
    #                                                '40A': '40A - Recurrent blood derived cancer – peripheral blood',
    #                                                '42A': '42A - Blood derived cancer- peripheral blood, post-treatment'})

    # df['Tissue Type'] = df['Tissue Code'].replace({'09A - Primary blood derived cancer - bone marrow': 'Bone Marrow',
    #                                                '03A - Primary blood derived cancer - peripheral blood': 'Peripheral Blood'})
    # df['Tumor Code'] = df['Tumor Code'].replace({'00': '00 - Non-cancerous tissue',
    #                                             '20': '20 - Acute myeloid leukemia (AML)',
    #                                              '21': '21 - Induction failure AML (AML-IF)'})
    return (df)


def clean_aml05(df, clinical_data_path='../Data/Raw_Data/Clinical_Data/'):
    """
    This module cleans and adjusts clinical data files from AML05 trial.

    Parameters:
    -----------
    df:object
        AML05 clinical outcome dataframe.
    
    clinical_data_path: str
        Path to the directory containing 'ClinicalData_AML05_JapaneseTrial_FromPaper.csv'.

    Returns:
    --------
        Cleaned and adjusted dataframe.

    """
    df = df.rename(columns={'gender': 'Gender', 'age': 'Age (years)',
                 'source': 'Sample Type'})[['Gender', 'Sample_Name',
                                            'Age (years)', 'Sample Type']]
    
    df = df.replace({'Unknown': np.nan, 'YES': 'Yes', 'NO': 'No', 'male': 'Male',
                     'female': 'Female', 'bone marrow': 'Bone Marrow'})

    # Keep only the second-last word of Sample_Name
    df['Sample_Name'] = df['Sample_Name'].str.split(' ').str[-2]

    # Load AML05 Clinical Data from Paper
    aml05 = pd.read_csv(
        clinical_data_path + 'ClinicalData_AML05_JapaneseTrial_FromPaper.csv', index_col=0)

    # All 15 samples are positive for FLT3 ITD according to paper's table 1
    aml05['FLT3 ITD'] = 'Yes'

    # Join AML05 Clinical Data to df
    df = df.join(aml05.drop(
        columns=['Age, y', 'Sex', 'Cluster']), on='Sample_Name')

    # Rename columns to match other datasets
    df = df.rename(columns={'FLT3-ITD AR': 'FLT3/ITD allelic ratio',
                            'Blast, %': 'BM Leukemic blasts (%)', 'WBC, 109/L': 'WBC Count (10⁹/L)',
                            'Outcome': 'Vital Status', 'Relapse': 'relapse.evnt', 'Gender': 'Sex',
                            'Sample_Name': 'Patient_ID', 'Sample Type': 'Tissue Type'})

    df['os.evnt'] = df['Vital Status'].map({'Alive': 0, 'Dead': 1})
    df['Age group (years)'] = pd.cut(x=df['Age (years)'].astype(float),
                                     bins=[-np.inf, 10, np.inf],
                                     labels=['<10', '≥10'])

    df['Race or ethnic group'] = 'Asian'  # Need to verify this
    df['Hispanic or Latino ethnic group'] = 'Not Hispanic or Latino'# Need to verify this
    df['Clinical Trial'] = 'Japanese AML05'
    df['Sample Type'] = 'Diagnosis'

    return (df)


def clean_beataml(df):

    df = df.rename(columns={'LLS_SampleID': 'Patient_ID', 'tissue':'Tissue', 'disease_state':'Sample Type',
                            'CEBPA_Biallelic': 'CEBPA mutation', 'FAB_BlastMorphology':'FAB',
                            'SpecificDxAtAcquisition':'Diagnosis','VitalStatus':'Vital Status',
                            'WBC_Count':'WBC Count (10⁹/L)', 'WHO_Fusion':'Gene Fusion'})
    
    df = df.replace({'Unknown': np.nan, 'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes'})
    df['Clinical Trial'] = 'Beat AML Consortium'
    return df


def clean_amltcga(df):

    df = df.rename(columns={'case_submitter_id':'Patient_ID', 'Age': 'Age (years)',
                            '%BM Blast':'BM Leukemic blasts (%)', 'WBC':'WBC Count (10⁹/L)',
                            '%PB Blast': 'Peripheral blasts (%)', 'Cytogenetics': 'Karyotype',
                            'Gene Fusions by RNA-Seq':'Gene Fusion', 
                            'Cytogenetic Classification':'Primary Cytogenetic Code',
                            'RISK (Molecular)':'Risk Group'})
    
    df = df.replace({'Unknown': np.nan, 'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes',
                     'M':'Male','F':'Female'})

    # Add `AML with` to the beginning of the `Diagnosis` column for each value
    df['Diagnosis'] = 'AML with ' + df['Molecular Classification'].astype(str)
    df['Clinical Trial'] = 'TCGA AML'
    df['Sample Type'] = 'Diagnosis'

    return df


def clean_nordic_all(df):

    df = df.rename(columns={'title':'Patient_ID', 'source':'Tissue', 'disease_state':'Sample Type',
                            'Karyotyping at diagnosisc':'Karyotype'})
    
    df = df.replace({'Unknown': np.nan, 'NA': np.nan,'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes',
                     'M':'Male','F':'Female', 'ALL diagnosis':'Diagnosis', 'ALL relapse':'Relapse'})

    # Make `Diagnosis` column by concatenating `immunophenotype` and `subtype` columns
    df['Diagnosis'] = df['immunophenotype'] + ' ' + df['subtype'] 
    df['Clinical Trial'] = 'NOPHO ALL92-2000'

    return df


def clean_mds_taml(df):

    df = df.rename(columns={'tissue': 'Tissue','diagnosis': 'Diagnosis', 'time_afer_diagnosis':'Sample Type',
                            'treatment':'Treatment'})
    
    df = df.replace({'Unknown': np.nan, 'NA': np.nan,'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes',
                     'M':'Male','F':'Female', 'bone marrow aspirate': 'Bone Marrow', 'DX':'Diagnosis',
                     'CTR':'Control (Healthy Donor)','POST TPH': 'Post Transplant'})
    
    df['Clinical Trial'] = 'CETLAM SMD-09 (MDS-tAML)'


    return df


def clean_all_graal(df):

    df = df.rename(columns={'source': 'Sample Type', 'age': 'Age (years)','gender':'Gender',
                            'diagnosis':'Diagnosis'})
    
    df = df.replace({'Unknown': np.nan, 'NA': np.nan,'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes',
                     'M':'Male','F':'Female', 'Leukemic Bone Marrow': 'Likely Diagnosis'})
    df['Clinical Trial'] = 'French GRAALL 2003–2005'


    return df


def clean_target_all(df):
      
    df['Clinical Trial'] = 'TARGET ALL'
    df['Sample Type'] = 'Unknown (TARGET-ALL)'

    return df

def label_control_samples(df_methyl, df):
    """
    This function labels control samples from the AML0531 clinical trial (GSE124413) as 'Bone Marrow Normal'
    and combines them with the clinical trial samples.
    """
    a = df_methyl[df_methyl['Batch'].isin(['GSE124413'])]
    b = df[df.index.isin(a.index)]
    control_0531 = a[~a.index.isin(b.index)]
    control_0531['Sample Type'] = 'Bone Marrow Normal'
    df_ = pd.concat(
        [df, control_0531['Sample Type'].to_frame()], axis=0, join='outer')
    return df_