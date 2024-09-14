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

import joblib
def load_model_and_impute_by_median(df, model_path):
    # Load the imputer model
    loaded_imputer = joblib.load(model_path)
    
    # Get the features that the imputer was trained on
    imputer_features = loaded_imputer.feature_names_in_
    
    # Transpose the DataFrame and select only the features the imputer knows about
    df_aligned = df.T.reindex(columns=imputer_features)
    
    # Transform the data
    imputed_data = loaded_imputer.transform(df_aligned)
    
    # Convert back to a DataFrame and transpose to original shape
    imputed_df = pd.DataFrame(imputed_data, columns=imputer_features, index=df_aligned.index).T
    
    return imputed_df.round(3).astype('float32')