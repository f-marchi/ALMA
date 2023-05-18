"""
This module implements functions for unsupervised learning algorithms.

"""

import numpy as np
import pandas as pd
import pacmap

def run_pacmap(x_train, x_test, n_components=2):
    """
    Run PaCMAP on the training dataset apply learned parameters to the train and test.

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