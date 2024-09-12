#python3.10
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm.notebook import tqdm

# Constants
MOUNT = Path('/mnt/e/')
PACMAP_REFERENCE_PATH = Path('../data/pacmap_reference.bed')
NANOPORE_PROCESSED_PATH = MOUNT / 'nanopore_processed/bed'

# Column names
BED_COLUMNS = [
    "chrom", "start_position", "end_position", "modified_base_code", "score",
    "strand", "start_position2", "end_position2", "color", "Nvalid_cov",
    "fraction_modified", "Nmod", "Ncanonical", "Nother_mod", "Ndelete",
    "Nfail", "Ndiff", "Nnocall"
]

def read_pacmap_reference(file_path):
    df = pd.read_csv(
        file_path, 
        sep='\t', 
        names=['chrm', 'start', 'end', 'name', 'score', 'strand'],
        usecols=['chrm', 'start', 'name'],  # Only read necessary columns
        dtype={'chrm': str, 'start': int, 'name': str}  # Specify dtypes for faster reading
    )
    df['coordinate'] = df['chrm'] + ':' + df['start'].astype(str)
    return df.set_index('coordinate')

def read_sample_data(file_path):
    return pd.read_csv(
        file_path, 
        sep='\t', 
        names=BED_COLUMNS,
        usecols=['chrom', 'start_position', 'modified_base_code', 'fraction_modified'],  # Only read necessary columnsb
        dtype={'chrom': str, 'start_position': int, 'modified_base_code': str, 'fraction_modified': float}  # Specify dtypes for faster reading
    )

def process_sample(pacmap_ref, sample_df, sample_name):
    # Filter and create coordinate column in one step
    sample_filtered = sample_df[sample_df['modified_base_code'] == 'm']
    sample_filtered['coordinate'] = sample_filtered['chrom'] + ':' + sample_filtered['start_position'].astype(str)
    sample_filtered.set_index('coordinate', inplace=True)
    
    # Merge data
    merged = pacmap_ref[['name']].join(sample_filtered[['fraction_modified']], how='inner')
    
    # Calculate beta values and prepare final Series
    beta_values = (merged['fraction_modified'] / 100).round(3)
    beta_values.index = merged['name']
    beta_values.name = sample_name
    
    return beta_values

def process_directory(directory_path, pacmap_reference):
    results = []
    bed_files = sorted(directory_path.glob('*.bed'))  # Sort bed files to ensure consistent order
    
    for bed_file in tqdm(bed_files, desc="Processing samples", unit="sample"):
        sample_name = bed_file.stem  # Use filename without extension as sample name
        sample_df = read_sample_data(bed_file)
        
        # Process sample
        result = process_sample(pacmap_reference, sample_df, sample_name)
        results.append(result)
    
    # Concatenate all results and sort columns
    return pd.concat(results, axis=1).T.sort_index(axis=1)

def bed_to_pickle():
    # Read pacmap reference
    pacmap_reference = read_pacmap_reference(PACMAP_REFERENCE_PATH)
    
    # Process all samples in the directory
    result = process_directory(NANOPORE_PROCESSED_PATH, pacmap_reference)
    
    dataset_title = f'{result.shape[0]}samples_{result.shape[1]}cpgs_nanopore_bvalues.pkl'

    # Save result as pickle file
    output_path = NANOPORE_PROCESSED_PATH / dataset_title
    result.to_pickle(output_path)
    
    return result