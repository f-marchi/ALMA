from sklearn.impute import SimpleImputer
import joblib
import numpy as np
import pandas as pd

def impute_and_save_by_median(df, output_path):
    """
    Impute NaN values in the given DataFrame, save the imputer model, and return the imputed DataFrame.
    
    Args:
    df (pandas.DataFrame): Input DataFrame with shape (features, samples)
    output_path (str): Path to save the imputer model
    strategy (str): Imputation strategy ('median', 'mean', or 'most_frequent')
    
    Returns:
    pandas.DataFrame: Imputed DataFrame with the same shape as input
    """

    imputer = SimpleImputer(strategy='median', missing_values=np.nan)
    imputer.fit(df.T)
    imputed_data = imputer.transform(df.T)
    
    # Convert back to a DataFrame with original shape
    imputed_df = pd.DataFrame(imputed_data.T, columns=df.columns, index=df.index)
    
    # Save the imputer model
    joblib.dump(imputer, output_path)
    print("CpGs imputed by median and model saved at", output_path)
    
    return imputed_df

def load_model_and_impute_by_median(df, model_path):
    """
    Load a saved imputer model and use it to impute NaN values in the given DataFrame.
    
    Args:
    df (pandas.DataFrame): Input DataFrame with shape (features, samples)
    model_path (str): Path to the saved imputer model
    
    Returns:
    pandas.DataFrame: Imputed DataFrame with the same shape as input
    """
    # Load the saved imputer model
    loaded_imputer = joblib.load(model_path)
    
    # Transform the data
    imputed_data = loaded_imputer.transform(df.T)
    
    # Convert back to a DataFrame with original shape
    imputed_df = pd.DataFrame(imputed_data.T, columns=df.columns, index=df.index).round(3).astype('float32')
    
    return imputed_df