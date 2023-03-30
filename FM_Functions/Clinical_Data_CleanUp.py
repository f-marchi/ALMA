"""
This module cleans and adjusts clinical data files from St. Jude and COG AML trials.

"""

import pandas as pd
import numpy as np

__author__ = 'Francisco Marchi, Lamba Lab, University of Florida'
__email__ = 'flourenco@ufl.edu'


def combine_and_index_clinicaldata():
    """
    This function combines all clinical data files into one dataframe and indexes it by the sample ID.

    Note: If this function breaks, it is likely that one or more original clinical data files
    have been modified or removed from assigned directory: (../Data/Raw_Data/Clinical_Data/)

    Returns:
    --------
    labels: pd.DataFrame
        FOur dataframes containing all clinical data from all trials and indexed by the sample ID.
        Order: labels_cog, labels_aml02, labels_aml08, labels_aml05
    """

    labels_0531 = pd.read_csv('../Data/Raw_Data/Clinical_Data/MarchiF_ClinicalData_AML0531_03P1.csv',
                              index_col=67)
    labels_0531['Sample Type'] = 'Diagnosis'

    labels2_1 = pd.read_csv('../Data/Raw_Data/Clinical_Data/MarchiF_ClinicalData_AAML1031.csv',
                            index_col=67).rename(columns={'Timepoint': 'Sample Type'})
    labels2_2 = pd.read_excel('../Data/Raw_Data/Clinical_Data/COGAAML1031.v3.11.9.21.xlsx',
                              index_col=0)

    labels_1031 = labels2_1.join(labels2_2[['Treatment Arm', 'KRAS', 'Protocol risk group classification']],
                                 how='left', on='Patient_ID')

    labels3_1 = pd.read_csv(
        '../Data/Raw_Data/GDC_TARGET_AML_Methyl450k/clinical_info_gdc_2022-06-28.tsv', sep='\t')
    COG_clinicaldata_all = pd.read_excel(
        '../Data/Raw_Data/Clinical_Data/COG_raw_clinicaldata.xlsx', index_col=1)

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

    labels_aml05 = pd.read_pickle(
        '../Data/Raw_Data/AML05_JapaneseTrial_GSE133986/GSE133986_GPL21145_meta_data.pkl'
    ).set_index('Sample_ID').rename(
        columns={'gender': 'Gender', 'age': 'Age (years)',
                 'source': 'Sample Type'})[['Gender', 'Sample_Name',
                                            'Age (years)', 'Sample Type']]

    # AML02

    labels5_1 = pd.read_excel('../Data/Raw_Data/Clinical_Data/JLamba-AML02-data-2019-02-12_JKL.Francisco.xlsx',
                              index_col=94)  # main clinical data file
    labels5_2 = pd.read_excel('../Data/Raw_Data/Clinical_Data/AML02_Clinical_Data_Mastersheet.xlsx',
                              index_col=0)  # only used to retrieve pLSC6 values
    labels5_3 = pd.read_excel('../Data/Raw_Data/Clinical_Data/AML02methylation_Lamba_July25.2014Summary_FM_Cleaned.xlsx',
                              index_col=1)  # used to merged clinical data with methyl

    # Drop column "DONOTUSE: ACCESSION_NBR" from labels5_1

    labels5_1 = labels5_1.drop(columns=['DONOTUSE: ACCESSION_NBR'])

    def clean_aml02(labels5_1, labels5_2, labels5_3):

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

    labels_aml02 = clean_aml02(labels5_1, labels5_2, labels5_3)

    labels6_1 = pd.read_excel('../Data/Raw_Data/Clinical_Data/AML08-clinical-AEs-IDs-2022-06-01.xlsx',
                              sheet_name=[0, 2],
                              index_col=0)

    labels6_2 = pd.read_csv('../Data/Raw_Data/Clinical_Data/AML02.AML08.PC.clin.rand.merge.N400.csv',
                            index_col=1)
    # Merge cleaned clinical data files with methylation data file
    labels6_1[0] = labels6_1[0].join(labels6_2['AGE'], how='left')
    labels_aml08 = labels6_1[2]['Meth.Sample.ID'].to_frame().join(labels6_1[0],
                                                                  how='inner').reset_index(
    ).set_index('Meth.Sample.ID')
    labels_aml08.columns = labels_aml08.columns.str.title()

    labels_cog = pd.concat(
        [labels_0531, labels_1031, labels_gdc_target], axis=0, join='outer')

    return (labels_cog, labels_aml02, labels_aml08, labels_aml05)


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
    df = df.replace({'Unknown': np.nan, 'YES': 'Yes', 'NO': 'No', 'male': 'Male',
                     'female': 'Female', 'bone marrow': 'Bone Marrow'})

    # Keep only the second-last word of Sample_Name
    df['Sample_Name'] = df['Sample_Name'].str.split(' ').str[-2]

    # Load AML05 Clinical Data from Paper
    aml05 = pd.read_csv(
        '../Data/Raw_Data/Clinical_Data/ClinicalData_AML05_JapaneseTrial_FromPaper.csv', index_col=0)

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
    # Need to verify this
    df['Hispanic or Latino ethnic group'] = 'Not Hispanic or Latino'
    df['Clinical Trial'] = 'AML05'
    df['Sample Type'] = 'Diagnosis'

    return (df)
