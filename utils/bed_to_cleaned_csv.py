import pandas as pd
import argparse
import os

def process_sample(input_path, sample_name):
    # Read the reference data
    pacmap_reference = pd.read_csv('bed/pacmap_reference.bed', sep='\t', header=None, 
                                   names=['chrm', 'start', 'end', 'name', 'score', 'strand'])

    # Read the input data
    df = pd.read_csv(os.path.join(input_path, 'bed', f'{sample_name}.bed'), 
                     sep='\s+',
                     names=["chrom", "start_position","end_position","modified base code","score",
                            'strand' ,"start position","end position", "color", "Nvalid_cov",
                            "fraction modified", "Nmod", "Ncanonical", "Nother_mod", "Ndelete",
                            "Nfail", "Ndiff", "Nnocall"])

    # Create 'coordinate' column for merging
    df['coordinate'] = df['chrom'].astype(str) + ':' + df['start_position'].astype(str)

    df_filtered = df[df['modified base code'].isin(['m'])].set_index('coordinate')

    pacmap_reference['coordinate'] = pacmap_reference['chrm'].astype(str) + ':' + pacmap_reference['start'].astype(str)
    pacmap_reference = pacmap_reference.set_index('coordinate')

    df_merged = pacmap_reference[['name']].join(df_filtered, how='inner')

    # Transform the fraction modified into beta values
    df_merged.loc[:, sample_name] = (df_merged['fraction modified'] / 100).round(3)

    df_processed = df_merged[['name', sample_name]].set_index('name').T

    # sort columns
    df_processed = df_processed.sort_index(axis=1)
    output_path = os.path.join(input_path, 'pacmap', f'{sample_name}_pacmap_bvalues.csv')

    # Save the final DataFrame
    df_processed.to_csv(output_path)
    print(f"Processed data saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Process sample data and generate pacmap b-values.")
    parser.add_argument("input_path", help="Path to the input directory")
    parser.add_argument("sample_name", help="Name of the sample")
    
    args = parser.parse_args()
    
    process_sample(args.input_path, args.sample_name)

if __name__ == "__main__":
    main()