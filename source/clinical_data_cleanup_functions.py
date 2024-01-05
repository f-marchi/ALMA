"""
This module cleans and adjusts clinical data files.

"""

import pandas as pd
import numpy as np

__author__ = 'Francisco Marchi, Lamba Lab, University of Florida'
__email__ = 'flourenco@ufl.edu'

##############################################################################################################

# Set paths to clinical data files
clinical_data_path ='../../Data/Clinical_Data/'
input_path_450k = '/mnt/d/MethylScore/Raw_Data/Methyl_Array_450k/'
input_path_EPIC = '/mnt/d/MethylScore/Raw_Data/Methyl_Array_EPIC/'

##############################################################################################################
# Merge and index clinical data files
##############################################################################################################

# COG AAML1031
def merge_index_1031():

    filepath1 = clinical_data_path + '/TARGET/TARGET-AML/TARGET_AML_ClinicalData_AML1031_20221108.xlsx'
    filepath2 = input_path_EPIC + '/GSE190931/sample_sheet_meta_data.pkl'

    # Load clinical data files
    labels_1031 = pd.read_excel(filepath1)
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
    labels_1031 = labels_1031.join(meta, how='right').reset_index().set_index('Sample_ID')

    # Rename columns
    labels_1031 = labels_1031\
        .drop(columns=['Gene Fusion'])\
        .rename(columns={'fusion':'Gene Fusion','timepoint':'Sample Type'})

    return labels_1031

# COG TARGET-AML
def merge_index_0531():

    dir      =  clinical_data_path + '/TARGET/TARGET-AML/'
    filepath1 = dir + 'TARGET_AML_ClinicalData_Discovery_20221108.xlsx'
    filepath2 = dir + 'TARGET_AML_ClinicalData_Validation_20221108.xlsx'
    filepath3 = dir + 'TARGET_AML_ClinicalData_LowDepthRNAseq_20221108.xlsx'
    filepath4 = input_path_EPIC + '/GSE124413/sample_sheet_meta_data.pkl'
    filepath5 = input_path_EPIC + '/GSE124413/GSE124413_series_matrix.csv'
    filepath6 = dir + 'gdc_sample_sheet.2023-05-13.tsv'

    # Load all clinical data files for 0531
    labels_0531_1 = pd.read_excel (filepath1, index_col=0)
    labels_0531_2 = pd.read_excel (filepath2, index_col=0)
    labels_0531_3 = pd.read_excel (filepath3, index_col=0)
    labels_gdc    = pd.read_csv   (filepath6, sep='\t')
    meta          = pd.read_pickle(filepath4)
    meta_matrix   = pd.read_csv   (filepath5)

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

    def adjust_labels_gdc(labels_gdc):
        '''
        This function adjusts the labels_gdc dataframe to match the labels_0531 dataframe.
        '''

        # Select only TARGET-AML samples
        labels_gdc = labels_gdc[labels_gdc['Project ID'] == 'TARGET-AML']

        # Select only rows in `File Name` column that contain `Grn.idat`
        labels_gdc = labels_gdc[labels_gdc['File Name'].str.contains('Grn.idat')]

        # Remove the `_Grn.idat` extension from `File Name` column
        labels_gdc['File Name'] = labels_gdc['File Name'].str.replace('_Grn.idat', '')

        # Select only relevant columns
        labels_gdc = labels_gdc[['File Name', 'Case ID', 'Sample ID', 'Sample Type']]\
                                .set_index('Case ID')
        
        # Rename columns
        labels_gdc = labels_gdc.rename(columns={'File Name': 'Sample_ID', 'Sample ID': 'Sample_ID_2'})

        return labels_gdc

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

    # Remove duplicate samples by keeping the row with fewer NaNs
    labels_0531 = remove_the_duplicate_samples_with_more_nulls()

    # Set index to `TARGET USI` and join with `Gene Fusion.1` column
    labels_0531 = labels_0531.set_index('TARGET USI').join(labels_0531_3[['Gene Fusion.1']], how='outer')

    # extract the last two words of `TARGET USI` by splitting on `-` 
    labels_0531['Tumor Code'] = labels_0531.index.str.split('-').str[1]
    labels_0531['Patient_ID'] = labels_0531.index.str.split('-').str[2]

    # Adjust labels_gdc dataframe   
    labels_gdc = adjust_labels_gdc(labels_gdc)

    # Join labels_0531 and labels_gdc dataframes
    labels_gdc = labels_0531.join(labels_gdc, how='right')\
                            .reset_index()\
                            .set_index('Sample_ID')

    # Adjust and clean metadata
    meta_cleaned = clean_meta(meta_matrix, meta)

    # drop duplicates in `Patient_ID` column
    labels_0531 = labels_0531.drop_duplicates(subset='Patient_ID', keep='first')\
        .reset_index()\
        .set_index('Patient_ID')

    # join with meta_cleaned dataframe
    labels_0531 = labels_0531.join(meta_cleaned, how='right')

    # drop columns that are not needed
    labels_0531 = labels_0531.drop(columns=['age', 'sex', 'GSM_ID'])

    # Rename values in `Sample Type` column
    labels_0531['Sample Type'] = labels_0531['Sample Type'].replace({'group: tumor': 'Diagnosis',
                                                                    'group: normal': 'Bone Marrow Normal'})

    # Set index to `Sample_ID` to match methylation samples
    labels_0531 = labels_0531.reset_index().set_index('Sample_ID')

    # Concat labels_0531 and labels_gdc dataframes
    labels = pd.concat([labels_0531, labels_gdc], axis=0, join='outer')

    # Combine `Gene Fusion` and `Gene Fusion.1` columns with comma while ignoring NaNs/empty strings
    labels['Gene Fusion'] = labels[['Gene Fusion', 'Gene Fusion.1']]\
    .apply(lambda x: ','.join(filter(lambda i: i is not None and i==i, x)), axis=1)

    # Replace '' with np.nan
    labels['Gene Fusion'] = labels['Gene Fusion'].replace({'': np.nan})

    # Drop `Gene Fusion.1` column
    labels = labels.drop(columns=['Gene Fusion.1'])

    return labels

# AML05
def merge_index_aml05():

    filepath1 = clinical_data_path + 'Japanese_AML05/AML05_sample_sheet_meta_data.pkl'

    df = pd.read_pickle(filepath1).iloc[:,:-1].set_index('Sample_ID')
    return (df)

# AML02
def merge_index_aml02():

    filepath1 = clinical_data_path + '/StJude_Trials/JLamba-AML02-data-2019-02-12_JKL.Francisco.xlsx'
    filepath2 = clinical_data_path + '/StJude_Trials/AML02_Clinical_Data_Mastersheet.xlsx'
    filepath3 = clinical_data_path + '/StJude_Trials/AML02methylation_Lamba_July25.2014Summary_FM_Cleaned.xlsx'


    labels1 = pd.read_excel(filepath1, index_col=94).drop(columns=['DONOTUSE: ACCESSION_NBR'])  # main clinical data file
    labels2 = pd.read_excel(filepath2, index_col=0)  # only used to retrieve pLSC6 values
    labels3 = pd.read_excel(filepath3, index_col=1)  # used to merged clinical data with methyl
    
    # Join contents of two columns into one
    labels3['Sample'] = labels3['Sentrix Barcode'].astype(str) + '_' + labels3['Sample Section'].astype(str)
    labels4 = labels3[labels3['type'] == 'AML02'].reset_index()
    
    # Replace - with _ in Sample ID
    labels4['Sample ID'] = labels4['Sample ID'].str.replace('-', '_')
    
    # Set index to Sample ID
    labels4 = labels4.set_index('Sample ID')

    # Combine all clinical data files together
    labels5 = labels1.join(labels2[['pLSC6_gb', 'pLSC6_Score', 'DNMT3B']],
                                how='left', on='U133A.Dx.ID').join(labels4['Sample'],
                                    how='inner', on='JL.id.x').reset_index().set_index('Sample')

    # Remove samples that Dr. Lamba requested to be removed
    labels6 = labels5[labels5['Meth450K.Array.ID'].isin(
        ['AML02_119_P022', 'AML02_147_P053', 'AML02_197_P022',
            'AML02_23_P022', 'AML02_39_P022', 'AML02_46_P022',
            'AML02_13_P022']).astype(int) == 0]

    return (labels6)

# AML08
def merge_index_aml08():

    filepath1 = clinical_data_path + '/StJude_Trials/AML08-clinical-AEs-IDs-2022-06-01.xlsx'
    filepath2 = clinical_data_path + '/StJude_Trials/AML02.AML08.PC.clin.rand.merge.N400.csv'

    labels1 = pd.read_excel(filepath1, sheet_name=[0, 2], index_col=0)
    labels2 = pd.read_csv(filepath2, index_col=1)

    # Merge cleaned clinical data files with methylation data file
    labels1[0] = labels1[0].join(labels2['AGE'], how='left')
    labels_aml08 = labels1[2]['Meth.Sample.ID'].to_frame().join(labels1[0],how='inner').\
                                                            reset_index().set_index('Meth.Sample.ID')
    labels_aml08.columns = labels_aml08.columns.str.title()

    return (labels_aml08)

# BeatAML
def merge_index_beataml():

    filepath1 = clinical_data_path + 'BeatAML/BEAT_AML_Raw clinical data_702.Samples.Vizome.xlsx'
    filepath2 = clinical_data_path + 'BeatAML/Beat_AML_561 Unique Samples_Cbioportal_OSHU.xlsx'
    filepath3 = input_path_EPIC + 'GSE159907/sample_sheet_meta_data.pkl'

    # Load the clinical data
    labels1 = pd.read_excel(filepath1, index_col=3)
    labels2 = pd.read_excel(filepath2)
    meta = pd.read_pickle(filepath3).iloc[:,:-1]

    # Create a new column with only the content inside [] from column 'Sample_Name'
    meta['LLS_SampleID'] = meta['Sample_Name'].str.extract(r"\[(.*?)\]", expand=False)

    # Set the index to the new column
    meta = meta[['tissue','disease_state','LLS_SampleID','Sample_ID']].set_index('LLS_SampleID')

    # Extract the last word from the 'Sample ID' column by splitting on the underscore
    labels2['LLS_SampleID'] = labels2['Sample ID'].str.split('_').str[-1]
    labels2 = labels2.set_index('LLS_SampleID')

    # Join the two dataframes
    labels_beataml = labels1.join(labels2, how='outer')

    # Join the two dataframes
    labels_beataml = meta.join(labels_beataml, how='left').reset_index().set_index('Sample_ID')

    return labels_beataml

# AML-TCGA
def merge_index_amltcga():

    filepath1 = clinical_data_path +'/TCGA_LAML/SuppTable01_NEJM2013_TCGA_AML.Paper_Mutation data.xlsx'
    filepath2 = clinical_data_path +'/TCGA_LAML/TCGA_Clinical data_NEJM 2013_Downloaded_Cbioportal_07-07-2020_Final.xlsx'
    filepath3 = clinical_data_path + 'TCGA_LAML/gdc_sample_sheet.2023-05-13.tsv'

    # Load the clinical data
    labels1    = pd.read_excel(filepath1, index_col=0).iloc[1:,:]
    labels2    = pd.read_excel(filepath2, index_col=1)
    labels_gdc = pd.read_csv(filepath3, sep='\t')

    def adjust_labels_gdc(labels_gdc):
        '''
        This function adjusts the labels_gdc dataframe to match the labels_0531 dataframe.
        '''

        # Select only TCGA-AML samples
        labels_gdc = labels_gdc[labels_gdc['Project ID'] == 'TCGA-LAML']

        # Select only rows in `File Name` column that contain `Grn.idat`
        labels_gdc = labels_gdc[labels_gdc['File Name'].str.contains('Grn.idat')]

        # Remove the `_Grn.idat` extension from `File Name` column
        labels_gdc['File Name'] = labels_gdc['File Name'].str.replace('_Grn.idat', '')

        # Select only relevant columns
        labels_gdc = labels_gdc[['File Name', 'Case ID', 'Sample ID', 'Sample Type']]\
                                .set_index('Case ID')
        
        # Rename columns
        labels_gdc = labels_gdc.rename(columns={'File Name': 'Sample_ID', 'Sample ID': 'Sample_ID_2'})

        return labels_gdc

    labels_gdc = adjust_labels_gdc(labels_gdc)

    # Add `TCGA-AB-` prefix to `TCGA Patient ID` column in labels1
    labels1['TCGA Patient ID'] = 'TCGA-AB-' + labels1['TCGA Patient ID'].astype(int).astype(str)

    # Set `TCGA Patient ID` column as index
    labels1 = labels1.reset_index().set_index('TCGA Patient ID')

    # Join labels1 and labels2
    labels = labels2.join(labels1, rsuffix='_dup', how='inner')

    # Remove duplicate columns (i.e. columns with `_dup` suffix)
    labels = labels.loc[:,~labels.columns.str.contains('_dup')]\
                .drop(['Sample Type'], axis=1)

    # Join labels and labels_gdc
    labels_amltcga = labels_gdc.join(labels, how='inner')\
                        .reset_index()\
                        .set_index('Sample_ID')
    
    labels_amltcga = labels_amltcga.rename(columns={'index': 'Patient_ID'})

    return labels_amltcga

# Nordic ALL
def merge_index_nordic_all():

    filepath1 = input_path_450k + 'GSE49031/sample_sheet_meta_data.pkl'
    filepath2 = clinical_data_path + 'Nordic_ALL/PMID_25729447_Supp_Clinical_Data.xlsx'

    # Load meta data from GSE49031
    meta = pd.read_pickle(filepath1).iloc[:,:-1].set_index('Sample_ID')

    # split meta `title` column by the last word
    meta['title'] = meta['title'].str.split().str[-1]

    # Set index to `title`
    meta = meta.reset_index().set_index('title')

    # Load clinical data from paper
    paper = pd.read_excel(filepath2, index_col=0,header=2,
                sheet_name='Table S7- Verification summary')\
                [['Karyotyping at diagnosisc']] # diagnosisc is not a typo, it's a superscript c at the end

    # Join meta and paper
    meta = meta.join(paper)

    # Reset index to `Sample_ID`
    meta = meta.reset_index().set_index('Sample_ID')

    return meta

# MDS_tAML
def merge_index_mds_taml():

    filepath1 = input_path_450k + 'GSE152710/sample_sheet_meta_data.pkl'
    filepath2 = clinical_data_path + 'MDS_tAML/PMID33446256_clinical_data_extracted_from_paper_supplement.xlsx'

    # Load meta data from GSE152710
    meta = pd.read_pickle(filepath1)\
            .iloc[:,:-1].set_index('description')[['tissue','response','Sample_ID','time_afer_diagnosis']]
        
    # Read clinical data from supplemental info from paper
    df = pd.read_excel(filepath2, index_col=1)
    
    # Join meta and df
    meta = meta.join(df, how='left')

    # Set index to `Sample_ID`
    meta = meta.reset_index().set_index('Sample_ID')        

    return meta
    
# Tcell_ALL_GRAAL
def merge_index_all_graal():

    filepath1 = input_path_EPIC + 'GSE147667/sample_sheet_meta_data.pkl'
    
    meta = pd.read_pickle(filepath1).iloc[:,:-1].set_index('Sample_ID')

    return meta
        
# GDC_TARGET-ALL
def merge_index_target_all():
    
    filepath1 = clinical_data_path +'/TARGET/TARGET-ALL-P3/TARGET_ALL_ClinicalData_Phase_III_20221108.xlsx'
    filepath2 = clinical_data_path + '/TARGET/TARGET-ALL-P3/gdc_sample_sheet.2023-05-13.tsv'

    # Load the clinical data
    labels1    = pd.read_excel(filepath1, index_col=0)
    labels_gdc = pd.read_csv(filepath2, sep='\t')

    def adjust_labels_gdc(labels_gdc):
        '''
        This function adjusts the labels_gdc dataframe to match the labels_0531 dataframe.
        '''

        # Select only TCGA-AML samples
        labels_gdc = labels_gdc[labels_gdc['Project ID'] == 'TARGET-ALL-P3']

        # Select only rows in `File Name` column that contain `Grn.idat`
        labels_gdc = labels_gdc[labels_gdc['File Name'].str.contains('Grn.idat')]

        # Remove the `_Grn.idat` extension from `File Name` column
        labels_gdc['File Name'] = labels_gdc['File Name'].str.replace('_Grn.idat', '')

        # Select only relevant columns
        labels_gdc = labels_gdc[['File Name', 'Case ID', 'Sample ID', 'Sample Type']]\
                                .set_index('Case ID')
        
        # Rename columns
        labels_gdc = labels_gdc.rename(columns={'File Name': 'Sample_ID', 'Sample ID': 'Sample_ID_2'})

        return labels_gdc

    labels_gdc = adjust_labels_gdc(labels_gdc)

    # Join labels and labels_gdc
    labels_target_all = labels_gdc.join(labels1, how='left')\
                        .reset_index()\
                        .set_index('Sample_ID')

    return labels_target_all

##############################################################################################################
# Clean clinical data files
##############################################################################################################

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
                            'Fullkaryotypestring': 'Karyotype'})
    
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
        {'Unknown': np.nan, 'Unkown': np.nan, 'Inevaluable': np.nan, -9999: np.nan, -9996: np.nan})
    
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

    df['Risk Group'] = df['Risk group'].replace({10: 'Low Risk', 30: 'High Risk'})
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

    df['Dx at Acquisition'] = df['Comment']
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


def clean_aml05(df):
    """
    This module cleans and adjusts clinical data files from AML05 trial.

    Parameters:
    -----------
    df:object
        AML05 clinical outcome dataframe.

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
        clinical_data_path + 'Japanese_AML05/ClinicalData_AML05_JapaneseTrial_FromPaper.csv', index_col=0)

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

    df['Race or ethnic group'] = 'Asian'  # TODO need to verify this
    df['Hispanic or Latino ethnic group'] = 'Not Hispanic or Latino'# TODO Need to verify this
    df['Clinical Trial'] = 'Japanese AML05'
    df['Sample Type'] = 'Diagnosis'

    df['CEBPA mutation'] = np.nan
    df['NPM mutation'] = np.nan
    df['Karyotype'] = np.nan
    df['Dx at Acquisition'] = np.nan

    return (df)


def clean_beataml(df):

    df = df.rename(columns={'LLS_SampleID': 'Patient_ID', 'tissue':'Tissue', 'disease_state':'Sample Type',
                            'CEBPA_Biallelic': 'CEBPA mutation','VitalStatus':'Vital Status',
                            'WBC_Count':'WBC Count (10⁹/L)', 'WHO_Fusion':'Gene Fusion', 
                            'Age at Diagnosis': 'Age (years)', 'NPM1 Consensus Call': 'NPM mutation',
                            'SpecificDxAtAcquisition': 'Dx at Acquisition',})
    
    df = df.replace({'Unknown': np.nan, 'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes'})
    df['Clinical Trial'] = 'Beat AML Consortium'

    return df


def clean_amltcga(df):

    df = df.rename(columns={'case_submitter_id':'Patient_ID', 'Age': 'Age (years)',
                            '%BM Blast':'BM Leukemic blasts (%)', 'WBC':'WBC Count (10⁹/L)',
                            '%PB Blast': 'Peripheral blasts (%)', 'Cytogenetics': 'Karyotype',
                            'Gene Fusions by RNA-Seq':'Gene Fusion',
                            'Cytogenetic Classification':'Primary Cytogenetic Code'})
    
    df = df.replace({'Unknown': np.nan, 'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes',
                     'M':'Male','F':'Female'})
    
    df['Risk Group'] = df['RISK (Molecular)'].replace({'Good': 'Low Risk',
                                                       'Intermediate': 'Standard Risk',
                                                       'Poor': 'High Risk',
                                                       'N.D.': np.nan})

    df['Dx at Acquisition'] = 'AML with ' + df['Molecular Classification'].astype(str)

    # Create `Yes` or `No` columns for mutations of interest
    df['CEBPA mutation'] = df['CEBPA'].isna().replace({True: 'No', False: 'Yes'})
    df['NPM mutation'] = df['NPM1'].isna().replace({True: 'No', False: 'Yes'})


    df['Clinical Trial'] = 'TCGA AML'
    df['Sample Type'] = 'Diagnosis'

    return df


def clean_nordic_all(df):

    df = df.rename(columns={'title':'Patient_ID', 'source':'Tissue', 'disease_state':'Sample Type',
                            'Karyotyping at diagnosisc':'Karyotype'})
    
    df = df.replace({'Unknown': np.nan, 'NA': np.nan,'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes',
                     'M':'Male','F':'Female', 'ALL diagnosis':'Diagnosis', 'ALL relapse':'Relapse',
                     'normal bone marrow':'Bone Marrow Normal (NOPHO ALL92-2000)',
                     'normal peripheral blood':'Peripheral Blood Normal (NOPHO ALL92-2000)',
                     })

    # Make `Diagnosis` column by concatenating `immunophenotype` and `subtype` columns
    df['Dx at Acquisition'] = df['immunophenotype'] + ' ' + df['subtype']
    
    df['Clinical Trial'] = 'NOPHO ALL92-2000'

    df['CEBPA mutation'] = np.nan
    df['NPM mutation'] = np.nan
    df['Karyotype'] = np.nan
    df['Gene Fusion'] = np.nan
    


    return df


def clean_mds_taml(df):

    df = df.rename(columns={'description': 'Patient_ID','tissue': 'Tissue', 'time_afer_diagnosis': 'Sample Type'})
    
    df = df.replace({'Unknown': np.nan, 'NA': np.nan,'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes',
                     'M':'Male','F':'Female', 'bone marrow aspirate': 'Bone Marrow', 'DX':'Diagnosis',
                     'CTR':'Bone Marrow Normal','POST TPH': 'Post Transplant'})
    
    df['Clinical Trial'] = 'CETLAM SMD-09 (MDS-tAML)'

    df['CEBPA mutation'] = np.nan
    df['NPM mutation'] = np.nan
    df['Karyotype'] = np.nan
    df['Gene Fusion'] = np.nan
    df['Dx at Acquisition'] = 'MDS-related or secondary myeloid neoplasms'

    return df


def clean_all_graal(df):

    df = df.rename(columns={'age': 'Age (years)','gender':'Gender',
                            'diagnosis':'Dx at Acquisition', 'Sample_Name': 'Patient_ID',
                            'tissue':'Tissue'})
    
    df = df.replace({'Unknown': np.nan, 'NA': np.nan,'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes',
                     'M':'Male','F':'Female'})
    
    
    df['Sample Type'] = 'Diagnosis'
    df['Clinical Trial'] = 'French GRAALL 2003–2005'

    df['CEBPA mutation'] = np.nan
    df['NPM mutation'] = np.nan
    df['Karyotype'] = np.nan
    df['Gene Fusion'] = np.nan

    return df


def clean_target_all(df):

    df = df.rename(columns={'Case ID': 'Patient_ID'})
      
    df['Clinical Trial'] = 'TARGET ALL'

    df['Age (years)'] = df["Age at Diagnosis in Days"].astype(
                        float).apply(lambda x: x/365)
    df['efs.time'] = df['Event Free Survival Time in Days'].astype(
                        float).apply(lambda x: x/365)
    df['os.time'] = df['Overall Survival Time in Days'].astype(
                        float).apply(lambda x: x/365)
    df['os.evnt'] = df['Vital Status'].map({'Alive': 0, 'Dead': 1})
    df['efs.evnt'] = df['First Event'].isin(
                        ['Relapse', 'Induction Failure', 'Death', 
                        'Death without Remission']).astype(int)

    df['CEBPA mutation'] = np.nan
    df['NPM mutation'] = np.nan
    df['Gene Fusion'] = np.nan 

    df['Dx at Acquisition'] = 'MPAL with ' + df['WHO ALAL Classification']


    return df


##############################################################################################################
