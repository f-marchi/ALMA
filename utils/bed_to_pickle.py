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
    return pd.read_csv(
        file_path, 
        sep='\t', 
        usecols=['chrm', 'start', 'name'],  # Only read necessary columns
        names=['chrm', 'start', 'end', 'name', 'score', 'strand'], 
        dtype={'chrm': str, 'start': int, 'name': str}
    ).assign(coordinate=lambda df: df['chrm'] + ':' + df['start'].astype(str)).set_index('coordinate')

def read_sample_data(file_path):
    return pd.read_csv(
        file_path, 
        sep='\t', 
        usecols=['chrom', 'start_position', 'modified_base_code', 'fraction_modified'],  # Only read necessary columns
        names=BED_COLUMNS, 
        dtype={'chrom': str, 'start_position': int, 'modified_base_code': str, 'fraction_modified': float}
    )

def process_sample(pacmap_ref, sample_df, sample_name):
    # Filter directly while creating the coordinate column
    sample_df['coordinate'] = sample_df['chrom'] + ':' + sample_df['start_position'].astype(str)
    
    # Perform merging operation
    merged = pacmap_ref[['name']].merge(sample_df[['coordinate', 'fraction_modified']], left_index=True, right_on='coordinate', how='inner')
    
    # Calculate beta values and prepare final Series
    beta_values = (merged['fraction_modified'] / 100).round(3)
    beta_values.index = merged['name']
    beta_values.name = sample_name
    
    return beta_values

def process_directory(directory_path, pacmap_reference):
    results = []
    bed_files = sorted(directory_path.glob('*.bed'))  # Sort bed files for consistent order
    
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