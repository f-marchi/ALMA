#!/usr/bin/env python
# coding: utf-8

# # Dimension Reduction with PaCMAP

# ## Where the data at?

# In[7]:


input_path = '../Data/Processed_Data/'
output_path = '../Data/Processed_Data/'


# ## Load Datasets

# In[8]:


import pandas as pd

x = pd.read_pickle(input_path+'x.pkl')
y = pd.read_csv(input_path+'y.csv', index_col=0)

print(
    f' Dataset (df) contains {x.shape[1]} rows (mC sites) and {x.shape[0]} columns (samples).')


# ## Train-Test Split

# Here we will split the data into a training/discovery and testing/validation set.
# 
# We will use ```y_train``` to denote the training set, and ```y_test``` to denote the testing set. 

# In[40]:


# Split train and test by clinical trial
y_train = y[~y['Clinical Trial'].isin(['AML02', 'AML08'])]
y_test = y[y['Clinical Trial'].isin(['AML02', 'AML08'])]

# Select samples in x that are in y_train
x_train = x.loc[y_train.index]
x_test = x.loc[y_test.index]


print(
    f"Discovery dataset (train) contains {x_train.shape[1]} rows (mC sites) and {x_train.shape[0]} columns (samples)")
print(
    f"\n{y_train['Clinical Trial'].value_counts(dropna=False).to_string()}\n")
print(
    f"Validation dataset (test) contains {x_test.shape[1]} rows (mC sites) and {x_test.shape[0]} columns (samples).")
print(f"\n{y_test['Clinical Trial'].value_counts(dropna=False).to_string()}\n")


# ## Batch Correction with pyCombat

# - __pyCombat__: a Python tool for batch effects correction in high-throughput molecular data using empirical Bayes methods
# 
# - __Github__: [https://epigenelabs.github.io/pyComBat/](https://epigenelabs.github.io/pyComBat/)
# 
# - __Paper__: [bioRxiv](https://doi.org/10.1101/2020.03.17.995431)

# In[43]:


from combat.pycombat import pycombat

# Correct batch effects in the training dataset
x_train2 = pycombat(x_train.T, y_train['Batch']).T

print('Succesfully corrected batch effects in the training dataset.')


# ## Dimension Reduction with PaCMAP

# - __PaCMAP__: Large-scale Dimension Reduction Technique Preserving Both Global and Local Structure
# 
# - __Github__: [https://github.com/YingfanWang/PaCMAP](https://github.com/YingfanWang/PaCMAP)
# 
# - __Paper__: [Journal of Machine Learning Research](https://jmlr.org/papers/v22/20-1061.html)

# In[56]:


import pacmap


def run_pacmap(x_train, x_test, n_components=2):
    """
    Run PaCMAP on the training dataset and apply the learned parameters to the train and test datasets.

    Parameters
    ----------
    x_train : pandas.DataFrame
        Training dataset.
    x_test : pandas.DataFrame
        Test dataset.
    n_components : int, optional
        Number of components. The default is 2.

    Returns
    -------
    embedding : numpy.ndarray
        Embedding of the training dataset.
    embedding_test : numpy.ndarray
        Embedding of the test dataset.

    """

    # Initialize PaCMAP. Note: hyperparameter tuning has been performed.
    reducer = pacmap.PaCMAP(n_components=n_components, n_neighbors=15,
                            MN_ratio=0.4, FP_ratio=16.0, random_state=42,
                            lr=0.1, num_iters=5000)

    # Fit (estimate) parameters to the training dataset to learn the embedding
    embedding = reducer.fit_transform(x_train)

    # Transform (apply) parameters to the test dataset
    embedding_test = reducer.transform(x_test, basis=x_train.copy())

    return embedding, embedding_test


embedding, embedding_test = run_pacmap(x_train2, x_test)


# ```{dropdown} Understanding Machine Learning
# 
# You may have noticed that we called two methods in the PaCMAP class: ```fit``` and ```transform```:
# 
# - __fit__: Learn the parameters of a model _from_ a dataset.
# 
# - __transform__: Apply the learned parameters _to_ a dataset.
# ```

# ## Save Embedding

# In[59]:


# Transform df to pandas dataframe format
embedding = pd.DataFrame(embedding, index=x_train2.index,
                         columns=['PaCMAP 1', 'PaCMAP 2'])
embedding_test = pd.DataFrame(embedding_test, index=x_test.index,
                              columns=['PaCMAP 1', 'PaCMAP 2'])

# Save embeddings
embedding.to_pickle(output_path+'embedding.pkl')
embedding_test.to_pickle(output_path+'embedding_test.pkl')

print(
    f'Successfuly saved {embedding.shape[0]} x_train samples and {embedding_test.shape[0]} x_test samples to {output_path}')


# ## Watermark

# In[47]:


get_ipython().run_line_magic('load_ext', 'watermark')


# In[58]:


# produce a list of the loaded modules
get_ipython().run_line_magic('watermark', '-v -p numpy,pandas,sklearn,combat,pacmap')

