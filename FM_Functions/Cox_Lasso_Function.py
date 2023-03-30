"""
This module implements a Cox Proportional Hazard regression with regularization.

"""

import pandas as pd
import numpy as np
import warnings

warnings.simplefilter(action='ignore')


def train_coxph_lasso(df, event, time, train_x=None, loops=10):
    """
    Implements a Cox Proportional Hazard regression with regularization and loops.
    For more information, please check dependencies.

    Parameters:
    ----------

    df: object
        A dataframe from training (discovery) dataset containg outcome/output labels.
        Note: labels need to be: os.evnt, os.time, efs.evnt, efs.time.
    event : object
        Identifier of column containing event indicator.
    time : object
        Identifier of column containing time to event.
    train_x: object
        A dataframe containing your features/variables.
    loops: int, default= 10
        The number of times the model should run at different random seeds.


    Returns:
    --------
        A dataframe containing the raw calculated coefficients.

    """

    # Import libraries for Cox-Lasso regression

    from sksurv.util import Surv
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sklearn.model_selection import GridSearchCV, KFold
    from sklearn.pipeline import make_pipeline
    from tqdm import tqdm

    print(f'Running Cox-Lasso through {loops} loops:')

    # Define x and y
    Y = Surv.from_dataframe(event, time, df.sort_index())
    X = train_x.sort_index()

    # Define range of alpha penalties to deploy
    # Note: alpha here means the constant that scales the lasso (L1) penalty in the cost function of the algorithm
    estimated_alphas = CoxnetSurvivalAnalysis(l1_ratio=1,).fit(X, Y).alphas_

    # Fit Cox-Lasso by calculating the ideal alpha and feature coefficients through 10-fold
    # cross-validation and random seed set to z number of loops (according to "loops")
    coefs = pd.DataFrame(index=X.columns)
    b = np.arange(loops)

    for i in tqdm(b):
        cv = KFold(n_splits=10, shuffle=True, random_state=i)
        gcv = GridSearchCV(make_pipeline(CoxnetSurvivalAnalysis(l1_ratio=1.0)),
                           param_grid={"coxnetsurvivalanalysis__alphas": [
                               [v] for v in estimated_alphas]},
                           cv=cv, n_jobs=-1).fit(X, Y)
        coefs = coefs.join(pd.DataFrame(gcv.best_estimator_.named_steps["coxnetsurvivalanalysis"].coef_,
                                        index=X.columns, columns=['coef_' + str(i)]))
    print('Cox-Lasso Trained Successfuly!')
    print(
        f'Selected alpha value: {gcv.best_estimator_.named_steps["coxnetsurvivalanalysis"].alphas}')
    # Select non-zero features (those not brought to 0 by the L1 (lasso) penalty)
    coefs['nonzero_count'] = coefs.astype(bool).sum(axis=1)
    coefs['nonzero_freq'] = coefs['nonzero_count']/b[-1]
    coefs.sort_values(by=['nonzero_freq'], ascending=False, inplace=True)
    return (coefs)


def set_cutoff(coefs, threshold):
    """
    Set feature inclusion threshold (only features present in [threshold]% of loops will be taken).

    Parameters:
    ----------
    coefs: object
        A dataframe containing the raw calculated coefficients. 
    threshold: float, default= 0.8
        Feature inclusion threshold (only features present in [threshold]% of loops will be taken)

    Returns:
    --------
    A dataframe containing the mean, non-zero calculated coefficients.

    """
    # Set feature inclusion threshold (only features present in 80% of loops will be taken)
    coefs2 = coefs[coefs['nonzero_freq'] > threshold].replace(
        0, np.NaN).drop(columns=['nonzero_count', 'nonzero_freq'])

    # Get the mean values of selected coefficients, ignoring 0s
    coef_mean = coefs2.mean(axis=1, skipna=True).sort_values(ascending=False)
    return (coef_mean)


def generate_coxph_score(coef_mean, x, df, score_name, train_test='train',
                         cutoff_train=0.5):
    """Generates a dataframe with score/Cox PH predictions

    Parameters:
    ----------
    coef_mean: object
        List of mean coeficients from CoxPH fit.
        Note: this value has to be a pandas series.
    df: object
        Dataframe to add your results to.
    x: object
        A dataframe of variables/features that will be used to calculate the score.
    score_name: 
        Name of your score/prediction column.
    train_test: str | float, default = 'train' 
        Use 'train' or float between 0 and 1 for cutoff percentile.
        If test, use pre-determined number for binary threshold/cutoff.
    cutoff: float, default=0.5
        Cutoff only matters if 'train' is selected in 'train_test' arg.
        cutoff is within (0,1) and defines cutoff percentile for
        categorical score.

    Returns:
    --------
        A dataframe with calculated quantitative and categorical score/model predictions
        and the value of the cutoff on the continuous variable.

    """

    # Calculate score in test(validation) data
    b = np.arange(len(coef_mean))
    df2 = pd.DataFrame()
    coef_mean = coef_mean.squeeze()  # Transforms coefs into pd.Series if pd.DataFrame

    for i in b:
        df2['coef_' + str(i)] = x[coef_mean.index[i]
                                  ].apply(lambda x: x * coef_mean[i])

    df3 = df2.iloc[:, (df2.shape[1] - coef_mean.shape[0]):df2.shape[1]]

    # Save calculated score as a separate column of the original dataframe
    df[score_name] = df3.sum(axis=1)

    # Determine train
    if train_test == 'train':
        i = np.quantile(df[score_name], cutoff_train)
    else:
        i = train_test

    print(f'Continuous score cut at the value of {round(i,4)}')

    # Binarize score
    df[score_name + ' Categorical'] = pd.cut(df[score_name],
                                             bins=[-np.inf, i, np.inf],
                                             labels=['Low', 'High'])

    df[score_name + '_cat_bin'] = pd.cut(df[score_name],
                                         bins=[-np.inf, i, np.inf],
                                         labels=[0, 1])

    return (df, i)
