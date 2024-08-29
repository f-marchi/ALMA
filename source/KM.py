"""
This module implements custom functions for several commonly used data visualization techniques.

"""

# Import libraries

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def draw_kaplan_meier(df, model_name, save_plot=False, figsize=(8, 10), 
                      add_risk_counts=False, save_survival_table=False,
                      trialname=None, show_ci=False):
    """
    Returns a Kaplan-Meier plot with:

        1. Hazard Ratio
        2. P-value
        3. Risk counts
        4. Survival table

    Parameters:
    ----------
    scorename: str
        Name of your model (and a column in df).
    df: object
        A dataframe containing:
            1. Continuous predictions of Cox Regression model under "scorename".
            2. efs/os information in the format of "efs.time" and "efs.evnt".
    save_plot: bool, default=False
        Set to True if you wish to save the plot.It will be saved under "../Figures/ForestPlot/"
    trialname: str
        Name of your clinical trial or dataset.
    scorename: str
        Name of your model.

    Returns:
    --------
        A magnificent double kaplan-meier figure.

    """
    # Import libraries for Kaplan Meier
    from lifelines.plotting import add_at_risk_counts
    from lifelines import KaplanMeierFitter
    from lifelines import CoxPHFitter

    # Set up the matplotlib figure
    sns.set_theme(style='white')
    f, ax = plt.subplots(2, 1, sharex=True, figsize=figsize)

    # Define survival curve categories
    groups = df[model_name]
    ix = (groups == 'High')

    # Fit the Kaplan Meier to each category of groups (kmf1 and kmf2)
    def surv_curves(i, t, e):

        T = df[t]
        E = df[e]
        kmf1 = KaplanMeierFitter()
        kmf1.fit(T[~ix], E[~ix], label='Low ' + model_name + ', n=' +
                 str(len(df[df[model_name].isin(['Low'])])))
        ax = kmf1.plot_survival_function(
            ax=i, show_censors=True, ci_show=show_ci)

        kmf2 = KaplanMeierFitter()
        kmf2.fit(T[ix], E[ix], label='High ' + model_name + ', n=' +
                 str(len(df[df[model_name].isin(['High'])])))
        ax = kmf2.plot_survival_function(
            ax=i, show_censors=True, ci_show=show_ci)

        try:
            # Calculate Hazard Ratio (HZ) and p-value (p)
            X_CPH = df[[model_name + '_int', t, e]]
            cph = CoxPHFitter()
            HZ = cph.fit(X_CPH, t, event_col=e)
            hz = HZ.hazard_ratios_[0]
            p = HZ.summary['p'][0]
            ci_lower = HZ.summary['exp(coef) lower 95%'][0]
            ci_upper = HZ.summary['exp(coef) upper 95%'][0]

            # Annotate HZ, CI, and p
            i.annotate(f'HR: {hz:.4f} (95% CI: {ci_lower:.2f}, {ci_upper:.2f})\n p-value: {p:.4f}',
               xy=(9.75, 0.085), xycoords='data',
               ha='right', va='center', fontsize=11,
               bbox={'boxstyle': 'round', 'facecolor': 'none',
                     'edgecolor': 'lightgray'})
        except:
            pass


        # Add risk counts below the graph
        if add_risk_counts == True:
            add_at_risk_counts(kmf1, kmf2, ax=i)

        # Save Survival Function Table
        if save_survival_table == True:
            surv1 = kmf1.survival_function_.join(kmf1.confidence_interval_)
            surv2 = kmf2.survival_function_.join(kmf2.confidence_interval_)
            surv3 = surv1.join(surv2, how='outer')
            surv3.to_csv('../../Figures/Kaplan_Meiers/KM_' + t + '_SurvivalTable_' +
                         model_name + '_' + trialname + '_' + str(len(df)) + '.csv')

        i.set_ylim(0, 1)
        i.set_ylabel("est. probability of survival $\hat{S}(t)$")

    surv_curves(i=ax[0], t='efs.time', e='efs.evnt')
    surv_curves(i=ax[1], t='os.time at 5y', e='Vital Status at 5y')

    ax[0].set_title('Event-Free Survival', loc='left',
                    pad=10, fontweight='bold', fontsize=10)
    ax[1].set_title('Overall Survival', loc='left', pad=10, fontweight='bold', fontsize=10)
    # Define Plot Specs
    plt.subplots_adjust(wspace=0, hspace=0.2)
    plt.suptitle(model_name + " in " + trialname + ", n=" + str(len(df)),
                 fontsize=11, y=0.94, fontweight='bold')
    plt.xlim(0, 10)
    plt.xlabel("time $t$ (years)")

    # Save plot figure
    if save_plot == True:
        plt.savefig('../Figures/Kaplan_Meiers/' + model_name + '_' + trialname + '_' + str(len(df)) + '.png',
                    bbox_inches='tight', dpi=300)

    return (plt.show())