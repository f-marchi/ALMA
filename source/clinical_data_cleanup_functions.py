"""
This module cleans and adjusts clinical data files.

"""

import pandas as pd
import numpy as np

__author__ = 'Francisco Marchi, Lamba Lab, University of Florida'
__email__ = 'flourenco@ufl.edu'


clinical_data_path ='../../Data/Clinical_Data/'
input_path_450k = '/mnt/d/MethylScore/Raw_Data/Methyl_Array_450k/'
input_path_EPIC = '/mnt/d/MethylScore/Raw_Data/Methyl_Array_EPIC/'

##############################################################################################################
# Merge and index clinical data files
##############################################################################################################

# COG AAML1031
def merge_index_1031():

    filepath1 = '/TARGET/TARGET-AML/TARGET_AML_ClinicalData_AML1031_20221108.xlsx'
    filepath2 = '../../Data/Raw_Data/Methyl_Array_EPIC/GSE190931/sample_sheet_meta_data.pkl'

    # Load clinical data files
    labels_1031 = pd.read_excel(clinical_data_path + filepath1)
    meta        = pd.read_pickle(clinical_data_path +filepath2)

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
    filepath1 = 'TARGET_AML_ClinicalData_Discovery_20221108.xlsx'
    filepath2 = 'TARGET_AML_ClinicalData_Validation_20221108.xlsx'
    filepath3 = 'TARGET_AML_ClinicalData_LowDepthRNAseq_20221108.xlsx'
    filepath4 = '../Data/Raw_Data/Methyl_Array_EPIC/GSE124413/sample_sheet_meta_data.pkl'
    filepath5 = '../Data/Raw_Data/Methyl_Array_EPIC/GSE124413/GSE124413_series_matrix.csv'
    filepath6 = 'gdc_sample_sheet.2023-05-13.tsv'

    # Load all clinical data files for 0531
    labels_0531_1 = pd.read_excel (dir + filepath1, index_col=0)
    labels_0531_2 = pd.read_excel (dir + filepath2, index_col=0)
    labels_0531_3 = pd.read_excel (dir + filepath3, index_col=0)
    labels_gdc    = pd.read_csv   (dir + filepath6, sep='\t')
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
    df = pd.read_pickle(
    clinical_data_path + '/Japanese_AML05/AML05_sample_sheet_meta_data.pkl'
        ).iloc[:,:-1].set_index('Sample_ID')
    return (df)

# AML02
def merge_index_aml02():

        labels5_1 = pd.read_excel(clinical_data_path + '/StJude_Trials/JLamba-AML02-data-2019-02-12_JKL.Francisco.xlsx',
                              index_col=94).drop(columns=['DONOTUSE: ACCESSION_NBR'])  # main clinical data file
        labels5_2 = pd.read_excel(clinical_data_path + '/StJude_Trials/AML02_Clinical_Data_Mastersheet.xlsx',
                              index_col=0)  # only used to retrieve pLSC6 values
        labels5_3 = pd.read_excel(clinical_data_path + '/StJude_Trials/AML02methylation_Lamba_July25.2014Summary_FM_Cleaned.xlsx',
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
        labels6_1 = pd.read_excel(clinical_data_path + '/StJude_Trials/AML08-clinical-AEs-IDs-2022-06-01.xlsx',
                                sheet_name=[0, 2],
                                index_col=0)

        labels6_2 = pd.read_csv(clinical_data_path + '/StJude_Trials/AML02.AML08.PC.clin.rand.merge.N400.csv',
                                index_col=1)
        # Merge cleaned clinical data files with methylation data file
        labels6_1[0] = labels6_1[0].join(labels6_2['AGE'], how='left')
        labels_aml08 = labels6_1[2]['Meth.Sample.ID'].to_frame().join(labels6_1[0],
                                                                    how='inner').reset_index(
        ).set_index('Meth.Sample.ID')
        labels_aml08.columns = labels_aml08.columns.str.title()
        return (labels_aml08)

# BeatAML
def merge_index_beataml():

    filepath1 = clinical_data_path +'BeatAML/BEAT_AML_Raw clinical data_702.Samples.Vizome.xlsx'
    filepath2 = clinical_data_path +'BeatAML/Beat_AML_561 Unique Samples_Cbioportal_OSHU.xlsx'
    filepath3 = '../Data/Raw_Data/Methyl_Array_EPIC/GSE159907/sample_sheet_meta_data.pkl'

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

    # Load meta data from GSE152710
    meta = pd.read_pickle('../Data/Raw_Data/Methyl_Array_450k/GSE152710/sample_sheet_meta_data.pkl')\
            .iloc[:,:-1].set_index('description')[['tissue','response','Sample_ID','time_afer_diagnosis']]
        
    # Read clinical data from supplemental info from paper
    df = pd.read_excel('../Data/Raw_Data/Clinical_Data/MDS_tAML/PMID33446256_clinical_data_extracted_from_paper_supplement.xlsx',
                           index_col=1)
    
    # Join meta and df
    meta = meta.join(df, how='left')

    # Set index to `Sample_ID`
    meta = meta.reset_index().set_index('Sample_ID')        

    return meta
    
# Tcell_ALL_GRAAL
def merge_index_all_graal():
        meta = pd.read_pickle('../Data/Raw_Data/Methyl_Array_EPIC/GSE147667/sample_sheet_meta_data.pkl')\
            .iloc[:,:-1].set_index('Sample_ID')
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
    df['WHO ALL 2022 Diagnosis'] = np.nan

    # ELN 2022 Diagnostic Annotation

    def classify_annotated_diagnosis_aml02(gene_fusion):
        mapping = {
        'Acute Megakaryoblastic Leukemia': 'AML with other rare recurring translocations', # adding FAB M7 AMKL to rare translocations assuming its either t(1;22) or CBFA2T3-GLIS2
        }
        
        for key, value in mapping.items():
            if key in gene_fusion:
                return value
            
    def classify_cytogenetic_code_aml02(gene_fusion2):
        mapping = {
        't (9;11)':'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        '11q23':'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        't (8;21)': 'AML with t(8;21)(q22;q22.1)/RUNX1::RUNX1T1',
        'inv (16)': 'AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22)/CBFB::MYH11',}
        
        for key, value in mapping.items():
            if key in gene_fusion2:
                return value
            
    def classify_karyotype_aml02(normal_samples):
        mapping = {
            'placeholder': 'placeholder'}
        
        for key, value in mapping.items():
            if key in normal_samples:
                return value
            

    def process_labels_aml02(df):
        
        df["ELN22_Diagnosis"] = (df["Primarydiagnosis"].astype(str).apply(classify_annotated_diagnosis_aml02))
        df['ELN22_Cytogenetics'] = df['Primary Cytogenetic Code'].astype(str).apply(classify_cytogenetic_code_aml02)
        df['ELN22_Karyotype'] = df['Karyotype'].astype(str).apply(classify_karyotype_aml02)

        df['ELN22 Combined Diagnoses'] = df[['ELN22_Diagnosis','ELN22_Cytogenetics','ELN22_Karyotype']]\
            .apply(lambda x: ','.join(filter(lambda i: i is not None and i==i, x)), axis=1)
            # Replace empty strings with NaN
        df['ELN22 Combined Diagnoses'] = df['ELN22 Combined Diagnoses'].replace('', np.nan)

        # Create `ELN 2022 Diagnosis` column by splitting `Combined Diagnosis` by comma and taking the first element
        df['ELN AML 2022 Diagnosis'] = df['ELN22 Combined Diagnoses'].str.split(',').str[0]

        # Drop columns created except for `ELN 2022 Diagnosis` and `Combined Diagnosis` columns
        df = df.drop(['ELN22_Karyotype','ELN22_Diagnosis','ELN22_Cytogenetics'], axis=1)
            
        return df

    df = process_labels_aml02(df)
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
    df['WHO ALL 2022 Diagnosis'] = np.nan

    # ELN 2022 and WHO 2022 Diagnostic Annotation

    def classify_annotated_diagnosis_aml08(gene_fusion):
        mapping = {
        'RUNX1-RUNX1T1': 'AML with t(8;21)(q22;q22.1)/RUNX1::RUNX1T1',
        't(8;21)':      'AML with t(8;21)(q22;q22.1)/RUNX1::RUNX1T1',
        'CBFB-MYH11':    'AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22)/CBFB::MYH11',
        'inv(16)':      'AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22)/CBFB::MYH11',
        't(16;16)':     'AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22)/CBFB::MYH11',
        'KMT2A':         'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'add(11)(q23)':  'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'MLL':           'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'PML-RARA':      'APL with t(15;17)(q24.1;q21.2)/PML::RARA',
        'DEK-NUP214':    'AML with t(6;9)(p23;q34.1)/DEK::NUP214',
        'MECOM':         'AML with inv(3)(q21.3q26.2) or t(3;3)(q21.3;q26.2)/MECOM-rearrangement',
        'NPM1-MLF1':     'AML with other rare recurring translocations',
        'PRDM16-RPN1':   'AML with other rare recurring translocations',
        'RBM15-MRTFA':   'AML with other rare recurring translocations',
        'RBM15-MKL1':    'AML with other rare recurring translocations',
        'NUP98':         'AML with other rare recurring translocations',
        'ETV6-MNX1':     'AML with other rare recurring translocations',
        'KAT6A-CREBBP':  'AML with other rare recurring translocations',
        'PICALM-MLLT10': 'AML with other rare recurring translocations',
        'FUS-ERG':       'AML with other rare recurring translocations',
        'RUNX1-CBFA2T3': 'AML with other rare recurring translocations',
        'CBFA2T3-GLIS2': 'AML with other rare recurring translocations',
        'BCR-ABL1':       'AML with t(9;22)(q34.1;q11.2)/BCR::ABL1'}
        
        for key, value in mapping.items():
            if key in gene_fusion:
                return value
            
    def classify_karyotype_aml08(normal_samples):
        mapping = {
            'placeholder': 'placeholder'}
        
        for key, value in mapping.items():
            if key in normal_samples:
                return value
            

    def process_labels_aml08(df):
        
        df["ELN22_Diagnosis"] = (df["Primarydiagnosis"].astype(str).apply(classify_annotated_diagnosis_aml08))
        df['ELN22_Karyotype'] = df['Karyotype'].astype(str).apply(classify_karyotype_aml08)

        df['ELN22 Combined Diagnoses'] = df[['ELN22_Diagnosis','ELN22_Karyotype']]\
            .apply(lambda x: ','.join(filter(lambda i: i is not None and i==i, x)), axis=1)
            # Replace empty strings with NaN
        df['ELN22 Combined Diagnoses'] = df['ELN22 Combined Diagnoses'].replace('', np.nan)

        # Create `ELN 2022 Diagnosis` column by splitting `Combined Diagnosis` by comma and taking the first element
        df['ELN AML 2022 Diagnosis'] = df['ELN22 Combined Diagnoses'].str.split(',').str[0]

        # Create `WHO 2022 Diagnosis` column by splitting `Combined Diagnosis` by comma and taking the first element
        df['WHO AML 2022 Diagnosis'] = df['ELN22 Combined Diagnoses'].str.split(',').str[0]

        # Drop columns created except for `ELN 2022 Diagnosis` and `Combined Diagnosis` columns
        df = df.drop(['ELN22_Karyotype','ELN22_Diagnosis'], axis=1)
            
        return df

    df = process_labels_aml08(df)

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



    # ELN 2022

    def classify_controls(normal_samples):
        mapping = {
            'Bone Marrow Normal': 'Otherwise-Normal Control',
            'Blood Derived Normal': 'Otherwise-Normal Control'}
        
        for key, value in mapping.items():
            if key in normal_samples:
                return value

    def classify_fusion(gene_fusion):
        mapping = {
        'RUNX1-RUNX1T1': 'AML with t(8;21)(q22;q22.1)/RUNX1::RUNX1T1',
        'CBFB-MYH11':    'AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22)/CBFB::MYH11',
        'KMT2A':         'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'add(11)(q23)':  'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'MLL':           'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'PML-RARA':      'APL with t(15;17)(q24.1;q21.2)/PML::RARA',
        'DEK-NUP214':    'AML with t(6;9)(p23;q34.1)/DEK::NUP214',
        'MECOM':         'AML with inv(3)(q21.3q26.2) or t(3;3)(q21.3;q26.2)/MECOM-rearrangement',
        'NPM1-MLF1':     'AML with other rare recurring translocations',
        'PRDM16-RPN1':   'AML with other rare recurring translocations',
        'RBM15-MRTFA':   'AML with other rare recurring translocations',
        'RBM15-MKL1':    'AML with other rare recurring translocations',
        'NUP98':         'AML with other rare recurring translocations',
        'ETV6-MNX1':     'AML with other rare recurring translocations',
        'KAT6A-CREBBP':  'AML with other rare recurring translocations',
        'PICALM-MLLT10': 'AML with other rare recurring translocations',
        'FUS-ERG':       'AML with other rare recurring translocations',
        'RUNX1-CBFA2T3': 'AML with other rare recurring translocations',
        'CBFA2T3-GLIS2': 'AML with other rare recurring translocations',
        'BCR-ABL1':       'AML with t(9;22)(q34.1;q11.2)/BCR::ABL1'}

        # Other uncharacterized abdnormalities present in the dataset but not in guidelines
        
        # 'CBFA2T3-GLIS3': 'AML with other rare recurring translocations',
        # 'PSIP1-NUP214':  'AML with other rare recurring translocations',
        # 'XPO1-TNRC18':   'AML with other rare recurring translocations', 
        # 'HNRNPH1-ERG':   'AML with other rare recurring translocations',
        # 'NIPBL-HOXB9':   'AML with other rare recurring translocations', 
        # 'SET-NUP214':    'AML with other rare recurring translocations', 
        # 'FLI1-IFIT2':    'AML with other rare recurring translocations', 
        # 'TCF4-ZEB2':     'AML with other rare recurring translocations',
        # 'MBTD1-ZMYND11': 'AML with other rare recurring translocations', 
        # 'FOSB-KLF6':     'AML with other rare recurring translocations', 
        # 'SFPQ-ZFP36L2':  'AML with other rare recurring translocations', 
        # 'RUNX1-LINC00478':'AML with other rare recurring translocations',
        # 'RUNX1-EVX1':     'AML with other rare recurring translocations',  
        # 'PSPC1-ZFP36L1':  'AML with other rare recurring translocations', 
        # 'EWSR1-FEV':      'AML with other rare recurring translocations',
        # 'STAG2-AFF2':     'AML with other rare recurring translocations', 
        # 'MYB-GATA1':      'AML with other rare recurring translocations', 
        # 'CBFA2T3-GLIS3':  'AML with other rare recurring translocations',
        # 'RUNX1-ZFPM2':    'AML with other rare recurring translocations', 
        # 'RUNX1-CBFA2T2':  'AML with other rare recurring translocations',
        # 'PIM3-BRD1':      'AML with other rare recurring translocations',
        # 'KAT6A-EP300':    'AML with other rare recurring translocations',
        # 'DOT1L-RPS15':    'AML with other rare recurring translocations',
        # 'FUS-FEV':        'AML with other rare recurring translocations',
        # 'KAT6A-NCOA2':    'AML with other rare recurring translocations',
        # 'JARID2-PTP4A1':  'AML with other rare recurring translocations',
        # 'FUS-FLI1':       'AML with other rare recurring translocations',     
        
        for key, value in mapping.items():
            if key in gene_fusion:
                return value

    def classify_cebpa(cebpa_mutation):
        mapping = {
            'Yes': 'AML with in-frame bZIP mutated CEBPA'}
        
        for key, value in mapping.items():
            if key in cebpa_mutation:
                return value

    def classify_npm(npm_mutation):
        mapping = {
            'Yes': 'AML with mutated NPM1',
        }

        for key, value in mapping.items():
            if key in npm_mutation:
                return value
            
    def classify_annotated_diagnosis(diagnosis):
        mapping = {
            'mutated NPM1': 'AML with mutated NPM1',
            'mutated CEBPA': 'AML with in-frame bZIP mutated CEBPA',
            'myelodysplasia-related changes': 'MDS-related or secondary myeloid neoplasms'
            }
        
        for key, value in mapping.items():
            if key in diagnosis:
                return value

    def process_labels_eln22(df):
        df['ELN 2022_Controls'] = df['Sample Type'].astype(str).apply(classify_controls)
        df['ELN 2022_Gene Fusion'] = df['Gene Fusion'].astype(str).apply(classify_fusion)
        df['ELN 2022_CEBPA'] = df['CEBPA mutation'].astype(str).apply(classify_cebpa)
        df['ELN 2022_NPM1'] = df['NPM mutation'].astype(str).apply(classify_npm)
        df['ELN 2022_Comment'] = df['Comment'].astype(str).apply(classify_annotated_diagnosis)

        df['ELN 2022 Combined Diagnoses'] = df[['ELN 2022_Controls','ELN 2022_Gene Fusion', 'ELN 2022_CEBPA', 'ELN 2022_NPM1', 'ELN 2022_Comment']]\
            .apply(lambda x: ','.join(filter(lambda i: i is not None and i==i, x)), axis=1)

        # Replace empty strings with NaN
        df['ELN 2022 Combined Diagnoses'] = df['ELN 2022 Combined Diagnoses'].replace('', np.nan)

        # Create `ELN 2022 Final Diagnosis` column by splitting `Combined Diagnosis` by comma and taking the first element
        df['ELN AML 2022 Diagnosis'] = df['ELN 2022 Combined Diagnoses'].str.split(',').str[0]

        # Drop columns created except for `ELN 2022 Final Diagnosis` and `Combined Diagnosis` columns
        df = df.drop(['ELN 2022_Controls','ELN 2022_Gene Fusion', 'ELN 2022_CEBPA', 'ELN 2022_NPM1', 'ELN 2022_Comment'], axis=1)
            
        return df

    # Process labels
    df = process_labels_eln22(df)

    # WHO 2022

    def classify_controls(normal_samples):
        mapping = {
            'Bone Marrow Normal': 'Otherwise-Normal Control',
            'Blood Derived Normal': 'Otherwise-Normal Control'}
        
        for key, value in mapping.items():
            if key in normal_samples:
                return value


    def classify_fusion(gene_fusion):
        mapping = {
        'RUNX1-RUNX1T1': 'AML with t(8;21)(q22;q22.1)/RUNX1::RUNX1T1',
        'CBFB-MYH11':    'AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22)/CBFB::MYH11',
        'KMT2A':         'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'add(11)(q23)':  'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'MLL':           'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'PML-RARA':      'APL with t(15;17)(q24.1;q21.2)/PML::RARA',
        'DEK-NUP214':    'AML with t(6;9)(p23;q34.1)/DEK::NUP214',
        'MECOM':         'AML with inv(3)(q21.3q26.2) or t(3;3)(q21.3;q26.2)/MECOM-rearrangement',
        'ETV6':          'AML with ETV6 fusion',
        'NPM1':          'AML with mutated NPM1',
        'RBM15-MKL1':    'AML with t(1;22)(p13.3;q13.1); RBM15::MKL1',
        'NUP98':         'AML with NUP98-fusion',
        'KAT6A-CREBBP':  'AML with t(8;16)(p11.2;p13.3); KAT6A::CREBBP',
        'FUS-ERG':       'AML with t(16;21)(p11;q22); FUS::ERG',
        'CBFA2T3-GLIS2': 'AML with CBFA2T3::GLIS2 (inv(16)(p13q24))',
        'BCR-ABL1':       'AML with t(9;22)(q34.1;q11.2)/BCR::ABL1'}

        # Other uncharacterized abdnormalities present in the dataset but not in guidelines
        # 'RUNX1-CBFA2T3': 'AML with other rare recurring translocations',
        # 'PRDM16-RPN1':   'AML with other rare recurring translocations',
        # 'PICALM-MLLT10': 'AML with other rare recurring translocations',
        # 'RBM15-MRTFA':   'AML with other rare recurring translocations',
        # 'CBFA2T3-GLIS3': 'AML with other rare recurring translocations',
        # 'PSIP1-NUP214':  'AML with other rare recurring translocations',
        # 'XPO1-TNRC18':   'AML with other rare recurring translocations', 
        # 'HNRNPH1-ERG':   'AML with other rare recurring translocations',
        # 'NIPBL-HOXB9':   'AML with other rare recurring translocations', 
        # 'SET-NUP214':    'AML with other rare recurring translocations', 
        # 'FLI1-IFIT2':    'AML with other rare recurring translocations', 
        # 'TCF4-ZEB2':     'AML with other rare recurring translocations',
        # 'MBTD1-ZMYND11': 'AML with other rare recurring translocations', 
        # 'FOSB-KLF6':     'AML with other rare recurring translocations', 
        # 'SFPQ-ZFP36L2':  'AML with other rare recurring translocations', 
        # 'RUNX1-LINC00478':'AML with other rare recurring translocations',
        # 'RUNX1-EVX1':     'AML with other rare recurring translocations',  
        # 'PSPC1-ZFP36L1':  'AML with other rare recurring translocations', 
        # 'EWSR1-FEV':      'AML with other rare recurring translocations',
        # 'STAG2-AFF2':     'AML with other rare recurring translocations', 
        # 'MYB-GATA1':      'AML with other rare recurring translocations', 
        # 'CBFA2T3-GLIS3':  'AML with other rare recurring translocations',
        # 'RUNX1-ZFPM2':    'AML with other rare recurring translocations', 
        # 'RUNX1-CBFA2T2':  'AML with other rare recurring translocations',
        # 'PIM3-BRD1':      'AML with other rare recurring translocations',
        # 'KAT6A-EP300':    'AML with other rare recurring translocations',
        # 'DOT1L-RPS15':    'AML with other rare recurring translocations',
        # 'FUS-FEV':        'AML with other rare recurring translocations',
        # 'KAT6A-NCOA2':    'AML with other rare recurring translocations',
        # 'JARID2-PTP4A1':  'AML with other rare recurring translocations',
        # 'FUS-FLI1':       'AML with other rare recurring translocations',     
        
        for key, value in mapping.items():
            if key in gene_fusion:
                return value

    def classify_cebpa(cebpa_mutation):
        mapping = {
            'Yes': 'AML with bZIP mutated CEBPA'}
        
        for key, value in mapping.items():
            if key in cebpa_mutation:
                return value

    def classify_npm(npm_mutation):
        mapping = {
            'Yes': 'AML with mutated NPM1',
        }

        for key, value in mapping.items():
            if key in npm_mutation:
                return value
            
    def classify_annotated_diagnosis(diagnosis):
        mapping = {
            'mutated NPM1': 'AML with mutated NPM1',
            'mutated CEBPA': 'AML with bZIP mutated CEBPA',
            'myelodysplasia-related changes': 'MDS-related or secondary myeloid neoplasms'
            }
        
        for key, value in mapping.items():
            if key in diagnosis:
                return value

    def process_labels_who21(df):
        df['WHO 2022_Controls'] = df['Sample Type'].astype(str).apply(classify_controls)
        df['WHO 2022_Gene Fusion'] = df['Gene Fusion'].astype(str).apply(classify_fusion)
        df['WHO 2022_CEBPA'] = df['CEBPA mutation'].astype(str).apply(classify_cebpa)
        df['WHO 2022_NPM1'] = df['NPM mutation'].astype(str).apply(classify_npm)
        df['WHO 2022_Comment'] = df['Comment'].astype(str).apply(classify_annotated_diagnosis)

        df['WHO 2022 Combined Diagnoses'] = df[['WHO 2022_Controls','WHO 2022_Gene Fusion', 'WHO 2022_CEBPA', 'WHO 2022_NPM1', 'WHO 2022_Comment']]\
            .apply(lambda x: ','.join(filter(lambda i: i is not None and i==i, x)), axis=1)

        # Replace empty strings with NaN
        df['WHO 2022 Combined Diagnoses'] = df['WHO 2022 Combined Diagnoses'].replace('', np.nan)

        # Create `WHO 2022 Final Diagnosis` column by splitting `Combined Diagnosis` by comma and taking the first element
        df['WHO AML 2022 Diagnosis'] = df['WHO 2022 Combined Diagnoses'].str.split(',').str[0]

        # Drop columns created except for `WHO 2022 Final Diagnosis` and `Combined Diagnosis` columns
        df = df.drop(['WHO 2022_Controls','WHO 2022_Gene Fusion', 'WHO 2022_CEBPA', 'WHO 2022_NPM1', 'WHO 2022_Comment'], axis=1)
            
        return df

    # Process labels
    df = process_labels_who21(df)

    return (df)


def clean_aml05(df, clinical_data_path='../Data/Raw_Data/Clinical_Data/Japanese_AML05/'):
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

    # ELN 2022

    def classify_fusion_aml05_eln22(gene_fusion):
        mapping = {
        'KMT2A':         'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'NPM1':          'AML with mutated NPM1',
        'NUP98':         'AML with other rare recurring translocations'}
        for key, value in mapping.items():
            if key in gene_fusion:
                return value
            
    # Rename `Other genetic alterations` column to `Gene Fusion`
    df = df.rename(columns={'Other genetic alterations': 'Gene Fusion'})

    df['ELN AML 2022 Diagnosis'] = df['Gene Fusion'].astype(str).apply(classify_fusion_aml05_eln22)

    # WHO 2022

    def classify_fusion_aml05_who21(gene_fusion):
        mapping = {
        'KMT2A':         'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'NPM1':          'AML with mutated NPM1',
        'NUP98':         'AML with NUP98-fusion'}
        for key, value in mapping.items():
            if key in gene_fusion:
                return value

    df['WHO AML 2022 Diagnosis'] = df['Gene Fusion'].astype(str).apply(classify_fusion_aml05_who21)

    return (df)


def clean_beataml(df):

    df = df.rename(columns={'LLS_SampleID': 'Patient_ID', 'tissue':'Tissue', 'disease_state':'Sample Type',
                            'CEBPA_Biallelic': 'CEBPA mutation','VitalStatus':'Vital Status',
                            'WBC_Count':'WBC Count (10⁹/L)', 'WHO_Fusion':'Gene Fusion', 
                            'Age at Diagnosis': 'Age (years)'})
    
    df = df.replace({'Unknown': np.nan, 'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes'})
    df['Clinical Trial'] = 'Beat AML Consortium'

    # ELN 2022

    def classify_annotated_diagnosis_beataml_eln22(gene_fusion):
        mapping = {
            "AML with mutated NPM1"                                         : "AML with mutated NPM1",
            "AML with myelodysplasia-related changes"                       : "MDS-related or secondary myeloid neoplasms",
            "AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22); CBFB-MYH11" : "AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22)/CBFB::MYH11",
            "AML with mutated CEBPA"                                        : "AML with in-frame bZIP mutated CEBPA",
            "Therapy-related myeloid neoplasms"                             : "MDS-related or secondary myeloid neoplasms",
            "PML-RARA"                                                      : "APL with t(15;17)(q24.1;q21.2)/PML::RARA",
            "AML with t(9;11)(p22;q23); MLLT3-MLL"                          : "AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement",
            "AML with t(8;21)(q22;q22.1); RUNX1-RUNX1T1"                    : "AML with t(8;21)(q22;q22.1)/RUNX1::RUNX1T1",
            "AML with inv(3)(q21q26.2) or t(3;3)(q21;q26.2); RPN1-EVI1"     : "AML with inv(3)(q21.3q26.2) or t(3;3)(q21.3;q26.2)/MECOM-rearrangement",
            "Mixed phenotype acute leukaemia, T/myeloid"                    : "Mixed phenotype acute leukemia T/myeloid",
            "Myeloid leukaemia associated with Down syndrome"               : "Myeloid leukaemia associated with Down syndrome",
        }
        for key, value in mapping.items():
            if key in gene_fusion:
                return value

    df["ELN AML 2022 Diagnosis"] = (df["SpecificDxAtAcquisition"]
                                    .astype(str)
                                    .apply(classify_annotated_diagnosis_beataml_eln22))
    # WHO 2022

    def classify_annotated_diagnosis_beataml_who21(gene_fusion):
        mapping = {
            "AML with mutated NPM1"                                         : "AML with mutated NPM1",
            "AML with myelodysplasia-related changes"                       : "MDS-related or secondary myeloid neoplasms",
            "AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22); CBFB-MYH11" : "AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22)/CBFB::MYH11",
            "AML with mutated CEBPA"                                        : "AML with bZIP mutated CEBPA",
            "Therapy-related myeloid neoplasms"                             : "MDS-related or secondary myeloid neoplasms",
            "PML-RARA"                                                      : "APL with t(15;17)(q24.1;q21.2)/PML::RARA",
            "AML with t(9;11)(p22;q23); MLLT3-MLL"                          : "AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement",
            "AML with t(8;21)(q22;q22.1); RUNX1-RUNX1T1"                    : "AML with t(8;21)(q22;q22.1)/RUNX1::RUNX1T1",
            "AML with inv(3)(q21q26.2) or t(3;3)(q21;q26.2); RPN1-EVI1"     : "AML with inv(3)(q21.3q26.2) or t(3;3)(q21.3;q26.2)/MECOM-rearrangement",
            "Mixed phenotype acute leukaemia, T/myeloid"                    : "Mixed phenotype acute leukemia T/myeloid",
            "Myeloid leukaemia associated with Down syndrome"               : "Myeloid leukemia associated with Down syndrome",
        }
        for key, value in mapping.items():
            if key in gene_fusion:
                return value

    df["WHO AML 2022 Diagnosis"] = (
        df["SpecificDxAtAcquisition"]
        .astype(str)
        .apply(classify_annotated_diagnosis_beataml_who21)
    )

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

    # Add `AML with` to the beginning of the `Diagnosis` column for each value
    df['Diagnosis'] = 'AML with ' + df['Molecular Classification'].astype(str)
    df['Clinical Trial'] = 'TCGA AML'
    df['Sample Type'] = 'Diagnosis'

    # ELN 2022 Diagnostic Annotation

    def classify_annotated_diagnosis_amltcga_eln22(gene_fusion):
        mapping = {
        'PML-RARA':         'APL with t(15;17)(q24.1;q21.2)/PML::RARA',
        'CBFB-MYH11':       'AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22)/CBFB::MYH11',
        'RUNX1-RUNX1T1':    'AML with t(8;21)(q22;q22.1)/RUNX1::RUNX1T1',
        'MLL':              'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'BCR-ABL1':         'AML with t(9;22)(q34.1;q11.2)/BCR::ABL1',
        'NUP98':            'AML with other rare recurring translocations'}
        for key, value in mapping.items():
            if key in gene_fusion:
                return value
            
    df['ELN AML 2022 Diagnosis'] = df['Molecular Classification']\
        .astype(str).apply(classify_annotated_diagnosis_amltcga_eln22)
    
    # WHO 2022 Diagnostic Annotation

    def classify_annotated_diagnosis_amltcga_who21(gene_fusion):
        mapping = {
        'PML-RARA':         'APL with t(15;17)(q24.1;q21.2)/PML::RARA',
        'CBFB-MYH11':       'AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22)/CBFB::MYH11',
        'RUNX1-RUNX1T1':    'AML with t(8;21)(q22;q22.1)/RUNX1::RUNX1T1',
        'MLL':              'AML with t(9;11)(p22;q23.3)/KMT2A-rearrangement',
        'NUP98':            'AML with NUP98-fusion'}
        for key, value in mapping.items():
            if key in gene_fusion:
                return value
            
    df['WHO AML 2022 Diagnosis'] = df['Molecular Classification']\
        .astype(str).apply(classify_annotated_diagnosis_amltcga_who21)
    

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
    df['Immunophenotype_Subtype'] = df['immunophenotype'] + ' ' + df['subtype'] 
    df['Clinical Trial'] = 'NOPHO ALL92-2000'

    # WHO 2022 Diagnostic Annotation

    def classify_controls(normal_samples):
        mapping = {
            'Bone Marrow Normal': 'Otherwise-Normal Control',
            'Peripheral Blood Normal': 'Otherwise-Normal Control'}
        
        for key, value in mapping.items():
            if key in normal_samples:
                return value


    def classify_fusion(immunophenotype_subtype):
        mapping = {
        'BCP-ALL HeH'           :'B-ALL with hyperdiploidy, high',
        'BCP-ALL t(12;21)'      :'B-ALL with t(12;21)(p13.2;q22.1); ETV6::RUNX1',         
        'BCP-ALL undefined'     :'B-ALL NOS',     
        'T-ALL T-ALL'           :'T-ALL NOS',
        'BCP-ALL non-recurrent' :'B-ALL NOS', 
        'BCP-ALL 11q23/MLL'     :'B-ALL with t(v;11q23.3); KMT2A-rearranged',       
        'BCP-ALL t(1;19)'       :'B-ALL with t(1;19)(q23;p13.3); TCF3::PBX1',       
        'BCP-ALL dic(9;20)'     :'B-ALL dic(9;20)',       
        'BCP-ALL t(9;22)'       :'B-ALL with t(9;22)(q34.1;q11.2); BCR::ABL1',       
        'BCP-ALL iAMP21'        :'B-ALL with iAMP21',        
        'BCP-ALL <45chr'        :'B-ALL with hypodiploidy',        
        'BCP-ALL >67chr'        :'B-ALL with hyperdiploidy'}

        # This criteria above needs to be revised if we choose to include ALL samples in the analysis 
        
        for key, value in mapping.items():
            if key in immunophenotype_subtype:
                return value
            
    def process_labels_nordic_all(df):
        df['WHO 2022_Controls'] = df['Sample Type'].astype(str).apply(classify_controls)
        df['WHO 2022_Immunophenotype_Subtype'] = df['Immunophenotype_Subtype'].astype(str).apply(classify_fusion)
        
        df['WHO 2022 Combined Diagnoses'] = df[['WHO 2022_Controls','WHO 2022_Immunophenotype_Subtype']]\
            .apply(lambda x: ','.join(filter(lambda i: i is not None and i==i, x)), axis=1)

        # Replace empty strings with NaN
        df['WHO 2022 Combined Diagnoses'] = df['WHO 2022 Combined Diagnoses'].replace('', np.nan)

        # Create `WHO 2022 Final Diagnosis` column by splitting `Combined Diagnosis` by comma and taking the first element
        df['WHO ALL 2022 Diagnosis'] = df['WHO 2022 Combined Diagnoses'].str.split(',').str[0]

        # Drop columns created except for `WHO 2022 Final Diagnosis` and `Combined Diagnosis` columns
        df = df.drop(['WHO 2022_Controls','WHO 2022_Immunophenotype_Subtype','WHO 2022 Combined Diagnoses'], axis=1)
            
        return df

    # Process labels
    df = process_labels_nordic_all(df)

    return df


def clean_mds_taml(df):

    df = df.rename(columns={'description': 'Patient_ID','tissue': 'Tissue', 'time_afer_diagnosis': 'Sample Type'})
    
    df = df.replace({'Unknown': np.nan, 'NA': np.nan,'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes',
                     'M':'Male','F':'Female', 'bone marrow aspirate': 'Bone Marrow', 'DX':'Diagnosis',
                     'CTR':'Control (Healthy Donor)','POST TPH': 'Post Transplant'})
    
    df['Clinical Trial'] = 'CETLAM SMD-09 (MDS-tAML)'

    # ELN 2022 and WHO 2022 Diagnostic Annotation

    def classify_controls_mds_taml(normal_samples):
        mapping = {
            'CTR': 'Otherwise-Normal Control'}
        
        for key, value in mapping.items():
            if key in normal_samples:
                return value

    def classify_annotated_diagnosis_mds_taml(gene_fusion):
        mapping = {
            "MDS": "MDS-related or secondary myeloid neoplasms",
            "AML": "MDS-related or secondary myeloid neoplasms",
        }
        for key, value in mapping.items():
            if key in gene_fusion:
                return value

    def process_labels_mds_taml(df):
        df['ELN22_Controls'] = df['Patient_ID'].astype(str).apply(classify_controls_mds_taml)
        df["ELN22_Diagnosis"] = (df["Cytological category group"].astype(str).apply(classify_annotated_diagnosis_mds_taml))

        df['ELN22 Combined Diagnoses'] = df[['ELN22_Controls','ELN22_Diagnosis']]\
            .apply(lambda x: ','.join(filter(lambda i: i is not None and i==i, x)), axis=1)
            # Replace empty strings with NaN
        df['ELN22 Combined Diagnoses'] = df['ELN22 Combined Diagnoses'].replace('', np.nan)

        # Create `ELN 2022 Diagnosis` column by splitting `Combined Diagnosis` by comma and taking the first element
        df['ELN AML 2022 Diagnosis'] = df['ELN22 Combined Diagnoses'].str.split(',').str[0]

        # Create `WHO 2022 Diagnosis` column by splitting `Combined Diagnosis` by comma and taking the first element
        df['WHO AML 2022 Diagnosis'] = df['ELN22 Combined Diagnoses'].str.split(',').str[0]

        # Drop columns created except for `ELN 2022 Diagnosis` and `Combined Diagnosis` columns
        df = df.drop(['ELN22_Controls','ELN22_Diagnosis'], axis=1)
            
        return df

    df = process_labels_mds_taml(df)

    return df


def clean_all_graal(df):

    df = df.rename(columns={'source': 'Tissue', 'age': 'Age (years)','gender':'Gender',
                            'diagnosis':'Diagnosis', 'Sample_Name': 'Patient_ID'})
    
    df = df.replace({'Unknown': np.nan, 'NA': np.nan,'YES': 'Yes', 'NO': 'No', 'n':'No', 'y':'Yes',
                     'M':'Male','F':'Female', 'Leukemic Bone Marrow': 'Bone Marrow',
                     'normal': 'Otherwise-Normal Control'})
    
    
    df['Sample Type'] = 'Diagnosis'
    df['Clinical Trial'] = 'French GRAALL 2003–2005'
        # WHO 2022 Diagnostic Annotation

    def classify_all_graal(diagnosis):
        mapping = {
        'Otherwise-Normal Control':'Otherwise-Normal Control',
        'T-ALL':'T-ALL NOS'}
        for key, value in mapping.items():
            if key in diagnosis:
                return value
            
    df['WHO ALL 2022 Diagnosis'] = df['Diagnosis'].astype(str).apply(classify_all_graal)

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
    # WHO 2022 Diagnostic Annotation

    def classify_controls_target_all(normal_samples):
        mapping = {
            'Bone Marrow Normal': 'Otherwise-Normal Control'}
        
        for key, value in mapping.items():
            if key in normal_samples:
                return value

    def classify_annotated_diagnosis_target_all(gene_fusion):
        mapping = {
            "T/M": "Mixed phenotype acute leukemia T/myeloid",
            'B/M': 'Mixed phenotype acute leukemia B/myeloid',
            'MLL': 'Mixed phenotype acute leukemia with t(v;11q23.3)/KMT2A-rearranged'
        }
        for key, value in mapping.items():
            if key in gene_fusion:
                return value

    def process_labels_target_all(df):
        df['WHO_Controls'] = df['Sample Type'].astype(str).apply(classify_controls_target_all)
        df["WHO_Diagnosis"] = (df["WHO ALAL Classification"].astype(str).apply(classify_annotated_diagnosis_target_all))

        df['WHO Combined Diagnoses'] = df[['WHO_Controls','WHO_Diagnosis']]\
            .apply(lambda x: ','.join(filter(lambda i: i is not None and i==i, x)), axis=1)
            # Replace empty strings with NaN
        df['WHO Combined Diagnoses'] = df['WHO Combined Diagnoses'].replace('', np.nan)

        # Create `WHO 2022 Diagnosis` column by splitting `Combined Diagnosis` by comma and taking the first element
        df['WHO ALL 2022 Diagnosis'] = df['WHO Combined Diagnoses'].str.split(',').str[0]

        # Drop columns created except for `WHO Final Diagnosis` and `Combined Diagnosis` columns
        df = df.drop(['WHO_Controls','WHO_Diagnosis','WHO Combined Diagnoses'], axis=1)
            
        return df

    df = process_labels_target_all(df)

    return df


##############################################################################################################
# Functions for merged dataset
##############################################################################################################

def process_df_labels(df):
    """
    Function to process a pandas dataframe, performing age categorization 
    and main disease classification.

    Args:
        df (pandas.DataFrame): Input dataframe.

    Returns:
        pandas.DataFrame: Processed dataframe with age categorized and main disease classified.
    """
    import numpy as np
    import pandas as pd

    def categorize_age(age):
        """
        Function to categorize age into a specific range.

        Args:
            age (int or float): The age to be categorized.

        Returns:
            str: A string representation of the age category.
        """
        if pd.isnull(age):
            return np.nan
        elif age < 5:
            return '0-5'
        elif age < 13:
            return '5-13'
        elif age < 39:
            return '13-39'
        elif age < 60:
            return '39-60'
        else:
            return '60+'

    def classify_main_disease(subtype):
        """
        Function to classify the main disease based on a given subtype.

        Args:
            subtype (str): The subtype of the disease.

        Returns:
            str: A string representation of the main disease.
        """
        mapping = {
            'AML': 'Acute myeloid leukemia (AML)',
            'ALL': 'Acute lymphoblastic leukemia (ALL)',
            'MDS': 'Myelodysplastic syndrome (MDS or MDS-like)',
            'Mixed phenotype acute leukemia': 'Mixed phenotype acute leukemia (MPAL)',
            'APL': 'Acute promyelocytic leukemia (APL)',
            'Otherwise-Normal Control': 'Otherwise-Normal (Control)'
        }

        for key, value in mapping.items():
            if key in subtype:
                return value

    def main_disease_class(df):
        """
        Function to classify the main disease and create a pathology class.

        Args:
            df (pandas.DataFrame): The dataframe to be processed.

        Returns:
            pandas.DataFrame: The processed dataframe with new columns for main disease and pathology class.
        """

        df['WHO_ALL'] = df['WHO ALL 2022 Diagnosis'].astype(str).apply(classify_main_disease)
        df['ELN_AML'] = df['ELN AML 2022 Diagnosis'].astype(str).apply(classify_main_disease)

        df['Hematopoietic Lineage'] = df[['ELN_AML', 'WHO_ALL']] \
            .apply(lambda x: ','.join(filter(lambda i: i is not None and i == i, x)), axis=1) \
            .replace('', np.nan)

        # Drop columns created except for `WHO Final Diagnosis` and `Combined Diagnosis` columns
        df = df.drop(['ELN_AML', 'WHO_ALL'], axis=1)

        return df

    # Convert 'Age (years)' to numeric, errors='coerce' will turn non-numeric data to NaN
    df['Age (years)'] = pd.to_numeric(df['Age (years)'], errors='coerce')

    # Then apply your function
    df['Age (group years)'] = df['Age (years)'].apply(categorize_age)
    
    # Process labels
    df = main_disease_class(df)
    
    # Create `WHO 2022 Diagnosis` column
    df['WHO 2022 Diagnosis'] = df[['WHO AML 2022 Diagnosis', 'WHO ALL 2022 Diagnosis']] \
            .apply(lambda x: ','.join(filter(lambda i: i is not None and i == i, x)), axis=1) \
            .replace('', np.nan)

    return df
