#!/usr/bin/env python
# coding: utf-8

# # Dimension Reduction with PaCMAP

# ## Where the data at?

# In[1]:


input_path = '../Data/Processed_Data/'
output_path = '../Data/Processed_Data/'


# ## Load Datasets

# In[2]:


import pandas as pd

x = pd.read_pickle(input_path+'x.pkl')
y = pd.read_csv(input_path+'y.csv', index_col=0)


# In[3]:


nanopore_sample = pd.read_pickle(input_path+'deci_flongle_runs.pkl')


# ## Test Train Split

# In[4]:


y['Clinical Trial'].value_counts(dropna=False)


# In[5]:


y_train = y[~y['Clinical Trial'].isin(['AML02','AML08'])]
y_test = y[y['Clinical Trial'].isin(['AML02','AML08'])]


# In[6]:


y_train['Clinical Trial'].value_counts(dropna=False)


# In[7]:


# Select samples in x that are in y_train
x_train = x.loc[y_train.index]
x_test = x.loc[y_test.index]


# In[8]:


x_train.shape, x_test.shape


# ## Batch Correction with pyCombat

# In[9]:


from combat.pycombat import pycombat
data_corrected = pycombat(x_train.T,y_train['Batch'])
x_train2 = data_corrected.T


# ## Dimension Reduction with PaCMAP

# ### Define

# In[10]:


import pacmap
reducer = pacmap.PaCMAP(n_components=2, n_neighbors=15, 
                            MN_ratio=0.4, FP_ratio=16.0,random_state=42,
                            lr=0.1, num_iters=5000)


# In[11]:


# import pacmap
# reducer = pacmap.PaCMAP(n_components=2, n_neighbors=15, 
#                             MN_ratio=0.5, FP_ratio=10.0,random_state=42,
#                             lr=0.1, num_iters=4500)


# In[12]:


embedding = reducer.fit_transform(x_train2)


# In[13]:


embedding_test = reducer.transform(x_test, basis=x_train2.copy())


# ## Save Embedding

# In[14]:


embedding = pd.DataFrame(embedding, index=x_train2.index, columns=['PaCMAP 1','PaCMAP 2'])
embedding.to_pickle(output_path+'embedding.pkl')


# In[15]:


embedding_test = pd.DataFrame(embedding_test, index=x_test.index, columns=['PaCMAP 1','PaCMAP 2'])
embedding_test.to_pickle(output_path+'embedding_test.pkl')


# ## Define Nanopore Sample

# In[16]:


kasumi2 = nanopore_sample.join(x_train.columns.to_frame(name='index'),how='right').set_index('index').T


# In[17]:


x_train_nano = pd.concat([x_train2,nanopore_sample['Deci0-Flongle-1'].dropna().to_frame(name='Deci0-Flongle-1').T],
                            axis=0,join='inner')


# In[18]:


reducer_nano = pacmap.PaCMAP(n_components=2, n_neighbors=15, 
                            MN_ratio=0.4, FP_ratio=16.0,random_state=42,
                            lr=0.1, num_iters=5000)
# Fit the training set
embedding2 = reducer_nano.fit_transform(x_train_nano)


# In[19]:


# embedding2_test2 = reducer_nano.transform(kasumi2.T['Deci0-Flongle-1'].dropna().to_frame().T, 
#                                             basis=x_train_nano.copy())


# In[20]:


embedding2 = pd.DataFrame(embedding2, index=x_train_nano.index, columns=['PaCMAP 1','PaCMAP 2'])
embedding2.to_pickle(output_path+'embedding_nano.pkl')

