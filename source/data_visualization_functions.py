"""
This module implements custom functions for several commonly used data visualization techniques.

"""

# Import libraries

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, accuracy_score



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
    surv_curves(i=ax[1], t='os.time', e='os.evnt')

    ax[0].set_title('Event-Free Survival', loc='left',
                    pad=10, fontweight='bold')
    ax[1].set_title('Overall Survival', loc='left', pad=10, fontweight='bold')
    # Define Plot Specs
    plt.subplots_adjust(wspace=0, hspace=0.2)
    plt.suptitle(model_name + " in " + trialname + ", n=" + str(len(df)),
                 fontsize='medium', y=0.94,
                 fontweight='bold')
    plt.xlim(0, 10)
    plt.xlabel("time $t$ (years)")

    # Save plot figure
    if save_plot == True:
        plt.savefig('../Figures/Kaplan_Meiers/' + model_name + '_' + trialname + '_' + str(len(df)) + '.png',
                    bbox_inches='tight', dpi=300)

    return (plt.show())


def draw_forest_plot(time, event, df, save_plot=False, trialname=None, model_name=None):

    """
    Generates a custom forest plot. The code in this function is not pretty but it just gets the job done.

    Parameters:
    ----------
    time: object
        List of mean coeficients from CoxPH fit.
        Note: this value has to be a pandas series.
    event: object
        Dataframe to add your results to.
    df: object
        A dataframe of variables/features that will be used to calculate the score.
    save_plot: bool, default=False
        Set to True if you wish to save the plot.It will be saved under "../Figures/ForestPlot/"
    trialname: str
        Name of your clinical trial or dataset.
    model_name: str
        Name of your model.

    Returns:
    --------
        A magnificent forest plot.

    """
    import myforestplot as mfp
    from tableone import TableOne
    import statsmodels.formula.api as smf
    import numpy as np
    
    fp = df[[model_name,
             'MRD 1 Status',
             'Risk Group',
             'FLT3 ITD',
             'Leucocyte counts (10⁹/L)',
             'Age group (years)',
             time, event]]

    event2 = event.replace('.', '_')
    time2 = time.replace('.', '_')

    if event[0] == 'o':
        event3 = 'OS'
    else:
        event3 = 'EFS'

    fp2 = fp.rename(columns={event: event2,
                             time: time2,
                             'MRD 1 Status': 'MRD_1_Status',
                             'FLT3 ITD': 'FLT3_ITD',
                             'Risk Group': 'Risk_Group',
                             'Leucocyte counts (10⁹/L)': 'WBC_count',
                             'Age group (years)': 'Age_group'})

    res = smf.phreg(formula=time2 + " ~ C("+model_name+",Treatment(reference='Low')) + C(MRD_1_Status) + C(Risk_Group,Treatment(reference='Low Risk')) + C(FLT3_ITD) + C(WBC_count) + C(Age_group)",
                    data=fp2, status=event2).fit()

    res2 = res.summary(xname=[model_name+'-High',
                              'MRD 1 Status-Positive',
                              'Risk Group-High Risk',
                              'Risk Group-Standard Risk',
                              'FLT3 ITD-Yes',
                              'Leucocyte counts (10⁹/L)-≥30',
                              'Age group (years)-≥10']).tables[1]

    res3 = res2.set_index(res2.index.str.split(pat='-', expand=True))

    mytable = TableOne(data=fp.drop(columns=[event, time]),
                       pval=False, missing=True, overall=True,
                       label_suffix=False, order={model_name: ['High'],
                                                  'MRD 1 Status': ['Positive'],
                                                  'Risk Group': ['High Risk', 'Standard Risk'],
                                                  'FLT3 ITD': ['Yes'],
                                                  'Leucocyte counts (10⁹/L)': ['≥30'],
                                                  'Age group (years)': ['≥10']}).tableone

    mytable2 = mytable.join(res3)

    mytable2["risk_pretty"] = mfp.add_pretty_risk_column(mytable2,
                                                         risk="HR",
                                                         lower='[0.025',
                                                         upper='0.975]',
                                                         fml=".2f"
                                                         )
    mytable3 = mytable2.reset_index(names=['category', 'item']).rename(columns={'HR': 'risk',
                                                                                '[0.025': 0,
                                                                                '0.975]': 1}).iloc[1:, :]

    mytable3['P>|t|'] = round(mytable3['P>|t|'], 4).replace(
        {np.nan: '', 0: '<0.0001'})

    plt.rcParams["font.size"] = 8
    fp = mfp.ForestPlot(df=mytable3,
                        ratio=[3, 3, 2],
                        fig_ax_index=[2],
                        dpi=300,
                        figsize=(9, 5),
                        vertical_align=True)
    fp.errorbar(index=2, errorbar_kwds=None)
    fp.axd[2].set_xlim([1, 8.5])
    fp.axd[2].set_xticks([0, 2, 4, 6, 8])
    fp.axd[2].set_xticklabels(labels=[0, 2, 4, 6, 8], fontdict={'fontsize': 8})
    fp.axd[2].set_xlabel("Hazard Ratio", fontsize=8)
    fp.axd[2].axvline(x=1, ymin=0, ymax=1.0, color="black", alpha=0.5)

    fp.axd[1].set_xlim([0.50, 1])
    fp.embed_cate_strings(1, "category", 0.5, header=trialname + " " + event3,
                          text_kwds=dict(fontweight="bold"),
                          header_kwds=dict(fontweight="bold"),
                          )
    fp.embed_strings(1, "item", 0.55, header="", replace={"age": ""})
    fp.embed_strings(1, "Overall", 0.86, header="n (%)")
    fp.embed_strings(3, "P>|t|", 0, header="P>|t|")
    fp.embed_strings(3, "risk_pretty", 0.4, header="Hazard Ratio (95% CI)")
    fp.horizontal_variable_separators()
    fp.draw_outer_marker(log_scale=False, scale=0.008, index=2)

    # Save plot figure
    if save_plot == True:
        plt.savefig('../Figures/Forest_Plots/' + model_name + '_' + trialname + '_' + str(len(df)) + '_' + event3 + '.png',
                    bbox_inches='tight', dpi=300)

    return (plt.show())


def draw_forest_plot_noMRD(time, event, df, save_plot=False, trialname=None, model_name=None):
    """
    Generates a custom forest plot. The code in this function is not pretty but it just gets the job done.

    Parameters:
    ----------
    time: object
        List of mean coeficients from CoxPH fit.
        Note: this value has to be a pandas series.
    event: object
        Dataframe to add your results to.
    df: object
        A dataframe of variables/features that will be used to calculate the score.
    save_plot: bool, default=False
        Set to True if you wish to save the plot.It will be saved under "../Figures/ForestPlot/"
    trialname: str
        Name of your clinical trial or dataset.
    model_name: str
        Name of your model.

    Returns:
    --------
        A magnificent forest plot.

    """
    import myforestplot as mfp
    from tableone import TableOne
    import statsmodels.formula.api as smf
    import numpy as np
    
    fp = df[[model_name,
             'Risk Group',
             'FLT3 ITD',
             'Leucocyte counts (10⁹/L)',
             'Age group (years)',
             time, event]]

    event2 = event.replace('.', '_')
    time2 = time.replace('.', '_')

    if event[0] == 'o':
        event3 = 'OS'
    else:
        event3 = 'EFS'

    fp2 = fp.rename(columns={event: event2,
                             time: time2,
                             'FLT3 ITD': 'FLT3_ITD',
                             'Risk Group': 'Risk_Group',
                             'Leucocyte counts (10⁹/L)': 'WBC_count',
                             'Age group (years)': 'Age_group'})

    res = smf.phreg(formula=time2 + " ~ C("+model_name+",Treatment(reference='Low')) + C(Risk_Group,Treatment(reference='Low Risk')) + C(FLT3_ITD) + C(WBC_count) + C(Age_group)",
                    data=fp2, status=event2).fit()

    res2 = res.summary(xname=[model_name+'-High',
                              'Risk Group-High Risk',
                              'Risk Group-Standard Risk',
                              'FLT3 ITD-Yes',
                              'Leucocyte counts (10⁹/L)-≥30',
                              'Age group (years)-≥10']).tables[1]

    res3 = res2.set_index(res2.index.str.split(pat='-', expand=True))

    mytable = TableOne(data=fp.drop(columns=[event, time]),
                       pval=False, missing=True, overall=True,
                       label_suffix=False, order={model_name: ['High'],
                                                  'Risk Group': ['High Risk', 'Standard Risk'],
                                                  'FLT3 ITD': ['Yes'],
                                                  'Leucocyte counts (10⁹/L)': ['≥30'],
                                                  'Age group (years)': ['≥10']}).tableone

    mytable2 = mytable.join(res3)

    mytable2["risk_pretty"] = mfp.add_pretty_risk_column(mytable2,
                                                         risk="HR",
                                                         lower='[0.025',
                                                         upper='0.975]',
                                                         fml=".2f"
                                                         )
    mytable3 = mytable2.reset_index(names=['category', 'item']).rename(columns={'HR': 'risk',
                                                                                '[0.025': 0,
                                                                                '0.975]': 1}).iloc[1:, :]

    mytable3['P>|t|'] = round(mytable3['P>|t|'], 4).replace(
        {np.nan: '', 0: '<0.0001'})

    plt.rcParams["font.size"] = 8
    fp = mfp.ForestPlot(df=mytable3,
                        ratio=[3, 3, 2],
                        fig_ax_index=[2],
                        dpi=300,
                        figsize=(9, 5),
                        vertical_align=True)
    fp.errorbar(index=2, errorbar_kwds=None)
    fp.axd[2].set_xlim([1, 8.5])
    fp.axd[2].set_xticks([0, 2, 4, 6, 8])
    fp.axd[2].set_xticklabels(labels=[0, 2, 4, 6, 8], fontdict={'fontsize': 8})
    fp.axd[2].set_xlabel("Hazard Ratio", fontsize=8)
    fp.axd[2].axvline(x=1, ymin=0, ymax=1.0, color="black", alpha=0.5)

    fp.axd[1].set_xlim([0.50, 1])
    fp.embed_cate_strings(1, "category", 0.5, header=trialname + " " + event3,
                          text_kwds=dict(fontweight="bold"),
                          header_kwds=dict(fontweight="bold"),
                          )
    fp.embed_strings(1, "item", 0.55, header="", replace={"age": ""})
    fp.embed_strings(1, "Overall", 0.86, header="n (%)")
    fp.embed_strings(3, "P>|t|", 0, header="P>|t|")
    fp.embed_strings(3, "risk_pretty", 0.4, header="Hazard Ratio (95% CI)")
    fp.horizontal_variable_separators()
    fp.draw_outer_marker(log_scale=False, scale=0.008, index=2)

    # Save plot figure
    if save_plot == True:
        plt.savefig('../Figures/Forest_Plots/' + model_name + '_' + trialname + '_' + str(len(df)) + '_' + event3 + '.png',
                    bbox_inches='tight', dpi=300)

    return (plt.show())


def draw_forest_plot_withBMblast(time, event, df, save_plot=False, trialname=None, model_name=None):
    """
    Generates a custom forest plot. The code in this function is not pretty but it just gets the job done.

    Parameters:
    ----------
    time: object
        List of mean coeficients from CoxPH fit.
        Note: this value has to be a pandas series.
    event: object
        Dataframe to add your results to.
    df: object
        A dataframe of variables/features that will be used to calculate the score.
    save_plot: bool, default=False
        Set to True if you wish to save the plot.It will be saved under "../Figures/ForestPlot/"
    trialname: str
        Name of your clinical trial or dataset.
    model_name: str
        Name of your model.

    Returns:
    --------
        A magnificent forest plot.

    """
    import myforestplot as mfp
    from tableone import TableOne
    import statsmodels.formula.api as smf
    import numpy as np
    
    fp = df[[model_name,
             'MRD 1 Status',
             'Risk Group',
             'FLT3 ITD',
             'Leucocyte counts (10⁹/L)',
             'BM Leukemic blasts (%) (grouped)', 
             'Age group (years)',
             time, event]]

    event2 = event.replace('.', '_')
    time2 = time.replace('.', '_')

    if event[0] == 'o':
        event3 = 'OS'
    else:
        event3 = 'EFS'

    fp2 = fp.rename(columns={event: event2,
                             time: time2,
                             'MRD 1 Status': 'MRD_1_Status',
                             'FLT3 ITD': 'FLT3_ITD',
                             'Risk Group': 'Risk_Group',
                             'Leucocyte counts (10⁹/L)': 'WBC_count',
                             'BM Leukemic blasts (%) (grouped)': 'BM_blasts',
                             'Age group (years)': 'Age_group'})

    res = smf.phreg(formula=time2 + " ~ C("+model_name+",Treatment(reference='Low')) + C(MRD_1_Status) + C(Risk_Group,Treatment(reference='Low Risk')) + C(FLT3_ITD) + C(WBC_count) + C(BM_blasts) + C(Age_group)",
                    data=fp2, status=event2).fit()

    res2 = res.summary(xname=[model_name+'-High',
                              'MRD 1 Status-Positive',
                              'Risk Group-High Risk',
                              'Risk Group-Standard Risk',
                              'FLT3 ITD-Yes',
                              'Leucocyte counts (10⁹/L)-≥30',
                              'BM Leukemic blasts (%) (grouped)->50',
                              'Age group (years)-≥10']).tables[1]

    res3 = res2.set_index(res2.index.str.split(pat='-', expand=True))

    mytable = TableOne(data=fp.drop(columns=[event, time]),
                       pval=False, missing=True, overall=True,
                       label_suffix=False, order={model_name: ['High'],
                                                  'MRD 1 Status': ['Positive'],
                                                  'Risk Group': ['High Risk', 'Standard Risk'],
                                                  'FLT3 ITD': ['Yes'],
                                                  'Leucocyte counts (10⁹/L)': ['≥30'],
                                                  'BM Leukemic blasts (%) (grouped)': ['>50'],
                                                  'Age group (years)': ['≥10']}).tableone

    mytable2 = mytable.join(res3)

    mytable2["risk_pretty"] = mfp.add_pretty_risk_column(mytable2,
                                                         risk="HR",
                                                         lower='[0.025',
                                                         upper='0.975]',
                                                         fml=".2f"
                                                         )
    mytable3 = mytable2.reset_index(names=['category', 'item']).rename(columns={'HR': 'risk',
                                                                                '[0.025': 0,
                                                                                '0.975]': 1}).iloc[1:, :]

    mytable3['P>|t|'] = round(mytable3['P>|t|'], 4).replace(
        {np.nan: '', 0: '<0.0001'})

    plt.rcParams["font.size"] = 8
    fp = mfp.ForestPlot(df=mytable3,
                        ratio=[3, 3, 2],
                        fig_ax_index=[2],
                        dpi=300,
                        figsize=(9, 5),
                        vertical_align=True)
    fp.errorbar(index=2, errorbar_kwds=None)
    fp.axd[2].set_xlim([1, 8.5])
    fp.axd[2].set_xticks([0, 2, 4, 6, 8])
    fp.axd[2].set_xticklabels(labels=[0, 2, 4, 6, 8], fontdict={'fontsize': 8})
    fp.axd[2].set_xlabel("Hazard Ratio", fontsize=8)
    fp.axd[2].axvline(x=1, ymin=0, ymax=1.0, color="black", alpha=0.5)

    fp.axd[1].set_xlim([0.50, 1])
    fp.embed_cate_strings(1, "category", 0.5, header=trialname + " " + event3,
                          text_kwds=dict(fontweight="bold"),
                          header_kwds=dict(fontweight="bold"),
                          )
    fp.embed_strings(1, "item", 0.55, header="", replace={"age": ""})
    fp.embed_strings(1, "Overall", 0.86, header="n (%)")
    fp.embed_strings(3, "P>|t|", 0, header="P>|t|")
    fp.embed_strings(3, "risk_pretty", 0.4, header="Hazard Ratio (95% CI)")
    fp.horizontal_variable_separators()
    fp.draw_outer_marker(log_scale=False, scale=0.008, index=2)

    # Save plot figure
    if save_plot == True:
        plt.savefig('../Figures/Forest_Plots/' + model_name + '_' + trialname + '_' + str(len(df)) + '_' + event3 + '.png',
                    bbox_inches='tight', dpi=300)

    return (plt.show())


def draw_boxplot(df, x, y, order, trialname, hue=None, save_plot=False, figsize=None):
    """
    Generates a custom box plot.

    Parameters:
    ----------
    df: object
        A dataframe containing the columns x, y,and hue.
    x: str
        Categorical variable for the x-axis of the plot.
    y: str
        Continuous variable for the y-axis of the plot.
    hue: str, default=None
        Optional variable to be used as hue.
    save_plot: bool, default=False
        Set to True if you wish to save the plot.It will be saved under "../Figures/Box_Plots/"
    figsize: tuple, default=None
        Tuple containing the figsize. Select None for automatic size selection.


    Returns:
    --------
        A magnificent box plot.

    """
    from statannotations.Annotator import Annotator

    # Set up the matplotlib figure
    sns.set_theme(style='white')
    plt.subplots(figsize=figsize)

    if order == 'auto':
        order2 = list(df[x].value_counts().index)
    else:
        order2 = order.copy()

    ax = sns.boxplot(y=y, x=x, data=df,
                     whis=[0, 100], width=.6, orient='v', order=order2, color='white')

    # Add in points to show each observation
    sns.stripplot(y=y, x=x, data=df, size=4, color=".3", linewidth=0, orient='v',
                  order=order2, hue=hue, palette='bright')

    # Tweak the visual presentation
    ax.xaxis.grid(False)

    sns.despine(trim=True, left=True)
    plt.title(y + ' by ' + x + ' in ' + trialname, fontsize='medium', y=1,
              fontweight='bold')
    plt.legend(loc='center right', bbox_to_anchor=(1.25, 0.5))

    if len(order2) <= 3:
        if len(order2) == 2:
            pairs = [(order2[0], order2[1])]
        elif len(order2) == 3:
            pairs = [(order2[0], order2[1]), (order2[0],
                                              order2[2]), (order2[1], order2[2])]

        # Annotate figure
        annotator = Annotator(ax, pairs, data=df, x=x, y=y, order=order2)
        annotator.configure(test='Kruskal', text_format='star', loc='inside',
                            comparisons_correction='bonf').apply_and_annotate()

        # Save plot figure
    if save_plot == True:
        plt.savefig('../Figures/Box_Plots/' + hue + '_' + trialname + '_' + str(len(df)) + '_' + x + '.png',
                    bbox_inches='tight', dpi=300)
    return (plt.show())


def draw_stacked_barplot(df, x, y, order, trialname, hue=None, save_plot=False, figsize=None, fontsize=10):
    """
    Generates a custom stacked bar plot.

    Parameters:
    ----------
    df: object
        A dataframe containing the columns x, y,and hue.
    x: str
        Categorical variable for the x-axis of the plot.
    y: str
        Continuous variable for the y-axis of the plot.
    hue: str, default=None
        Optional variable to be used as hue.
    save_plot: bool, default=False
        Set to True if you wish to save the plot. It will be saved under "../Figures/Bar_Plots/"
    figsize: tuple, default=None
        Tuple containing the figsize. Select None for automatic size selection.


    Returns:
    --------
        A magnificent bar plot.

    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import numpy as np
    import pandas as pd

    # Set up the matplotlib figure
    sns.set_theme(style='white')
    plt.subplots(figsize=figsize)

    if order == 'auto':
        order2 = list(df[x].value_counts().index)
    else:
        order2 = order.copy()

    # Count the occurrences of each hue within each x value
    hue_counts = df.groupby([x, hue])[y].count().unstack(fill_value=0)
    # Convert these counts to proportions
    hue_props = hue_counts.divide(hue_counts.sum(axis=1), axis=0)
    
    # Generate a stacked bar plot
    hue_props.loc[order2, :].plot(kind='bar', stacked=True, ax=plt.gca())

    # Annotate the bars with their percentage values
    for rect in plt.gca().patches:
        height = rect.get_height()
        plt.gca().text(rect.get_x() + rect.get_width() / 2, rect.get_y() + height / 2,
                       '{:.1f}%'.format(height * 100), ha='center', va='center', color='white', fontsize=fontsize)

    plt.gca().yaxis.set_major_formatter(mticker.PercentFormatter(1))

    # Tweak the visual presentation
    plt.gca().xaxis.grid(False)

    # Turn tick labels 90 degrees
    plt.xticks(rotation=0,ha='center')

    plt.title(y + ' by ' + x + ' in ' + trialname, fontsize='medium', y=1)
    plt.legend(loc='center right', bbox_to_anchor=(1.25, 0.5))

    # Save plot figure
    if save_plot == True:
        plt.savefig('../Figures/Bar_Plots/' + '_'.join(hue_props.columns.tolist()) + '_' + trialname + '_' + str(len(df)) + '_' + x + '.png',
                    bbox_inches='tight', dpi=300)

    return plt.show()


def draw_heatmaps(fig_title, t1, t2, df1, df2, fig_number, save_plot=False, figsize=(9, 7)):
    """
    Generates a double-heatmap based on two datasets.

    Parameters:
    ----------
    fig_title: str
        The title of your figure.
    t1: str
        The title of your first panel.
    t2: str
        The title of your second panel.
    df1: object
        A dataframe containing the data for panel 1.
    df2: object
        A dataframe containing the data for panel 2.
    fig_number: int
        The number of your figure.
    save_plot: bool, default=False
        Set to True if you wish to save the plot.
        Note: It will be saved under "../Figures/Heatmaps/"
    figsize: tuple, default=(9, 7)
        Tuple containing the figsize. Select None for automatic size selection.


    Returns:
    --------
        A magnificent double-heatmap.

    """

    _, (ax1, ax2, axcb) = plt.subplots(1, 3,
                                       gridspec_kw={
                                           'width_ratios': [1, 1, 0.06]},
                                       figsize=figsize)

    def draw_ax(df, ax, y_label, t):

        sns.heatmap(df, ax=ax, xticklabels=False, yticklabels=y_label,
                    vmax=1, vmin=0, center=0.5, cbar_ax=axcb)
        ax.set_xlabel('Patient samples, n='+str(df.shape[1]), fontsize=10)
        ax.set_title(t, fontsize=11, loc='center')
        return ()

    draw_ax(df1, ax=ax1, y_label=True, t=t1)
    draw_ax(df2, ax=ax2, y_label=False, t=t2)

    plt.suptitle(t='Fig ' + str(fig_number) +
                 '. ' + fig_title,
                 fontsize='medium', fontweight='bold',
                 ha='right', y=0.98)

    # Adjusts fontsize for y labels
    ax1.set_yticklabels(ax1.get_yticklabels(), fontsize=10)
    axcb.set_yticklabels(axcb.get_yticklabels(), fontsize=8)

    if save_plot == True:
        plt.savefig('../Figures/Heatmaps/' + fig_title,
                    bbox_inches='tight', dpi=300)

    return (plt.show())


def draw_scatterplot(df_train, df_test, x, y, hue, s=25, save_plot=False, figsize=(7, 7)):
    """
    Generates a custom scatterplot.

    Parameters:
    ----------
    df_train: object
        A dataframe containing the columns x, y, and hue with training/discovery data.
    df_test: object
        A dataframe containing the columns x, y, and hue with testing/validation data.
    x: str
        Continuous variable for the x-axis of the plot.
        Note: NaN values are not counted.
    y: str
        Continuous variable for the y-axis of the plot.
    hue: str, default=None
        Optional variable to be used as hue.
    s: int, default=25
        Determines the diameter of the data points.
    save_plot: bool, default=False
        Set to True if you wish to save the plot.
        Note: It will be saved under "../Figures/Box_Plots/"
    figsize: tuple, default=(7,7)
        Tuple containing the figsize. Select None for automatic size selection.


    Returns:
    --------
        A magnificent scatterplot.

    """
    import scipy.stats as stats

    f, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=figsize)
    sns.despine(f, left=False, bottom=False)

    def draw_ax(df, ax, annot, legend):
        df = df[df[y] > -np.inf]
        sns.scatterplot(x=x, y=y, hue=hue, s=s, sizes=(1, 8), linewidth=0, alpha=0.9,
                        data=df, ax=axs[ax], legend=legend)

        r, p = stats.pearsonr(df[x], df[y])

        plt.annotate(fr'$\rho$: {r:.5f}' +
                     f'\np-value: {p:.4f}',
                     xy=annot, xycoords='figure fraction',
                     ha='right', va='center',
                     bbox={'boxstyle': 'round', 'facecolor': 'none',
                           'edgecolor': 'lightgray'})
        return ()

    draw_ax(df_train, ax=0, annot=(0.85, 0.56), legend='auto')
    draw_ax(df_test, ax=1, annot=(0.85, 0.14), legend=False)

    axs[0].set_title(' Discovery (AML02 trial), n=' +
                     str(len(df_train)), loc='center', pad=5, fontsize=11)
    axs[1].set_title(' Validation (COG trials), n=' +
                     str(len(df_test)), loc='center', pad=5, fontsize=11)

    plt.xlabel('Fraction of ' + x)
    plt.ylabel(y)
    plt.suptitle(y + ' by ' + x, fontsize='medium', y=0.95,
                 fontweight='bold')

    if save_plot == True:
        plt.savefig('../Figures/Scatterplots/' +
                    y+'_'+x+'.png',
                    bbox_inches='tight', dpi=300)

    return (plt.show())


def draw_pacmap(score, labels, hue=None, test_sample=None, panel_size=50, s=10, legend=True):
    """
    Generates a custom scatterplot.

    Parameters:
    ----------
    score: object
        A dataframe containing the columns x, y, and hue with training/discovery data.
    labels: object
        A dataframe containing the columns x, y, and hue with testing/validation data.
    hue: str, default=None
            Optional variable to be used as hue.
    panel_size: int, default=40
        Determines the size of the panel.
    test_sample: tuple, default=None
        Optional tuple containing the coordinates of the test sample.
        Note: The tuple must be in the form (x, y).
    s: int, default=10
        Determines the diameter of the data points.
    legend: bool, default=True
        Set to False if you wish to hide the legend.




    Returns:
    --------
    A magnificent scatterplot.
                """

    sns.set_theme(style="white")

    # Define variables
    score2 = score[:, 0:2]
    xs = score2[:, 0]
    ys = score2[:, 1]

    # Define scatterplot
    plt.subplots(figsize=(7, 5))

    sns.scatterplot(data=labels, x=xs, y=ys,
                    s=s, hue=hue,
                    linewidth=0, alpha=0.8, legend=legend)
    if test_sample != None:
        plt.scatter(test_sample[0], test_sample[1],
                    s=100, c='black', marker='*')

    # Define plot specs

    plt.xlabel("PaCMAP 1")
    plt.ylabel("PaCMAP 2")
    plt.tight_layout()
    plt.xlim(-panel_size, panel_size)
    plt.ylim(-panel_size, panel_size)
    # plt.grid(True)
    # Put the legend out of the figure
    if legend == True:
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # Define plot specs

    if hue != None:
        plt.title("AML Methylome Map by " + hue + ", n= " + str(len(score2)),
                  fontsize=12)
        plt.savefig('../Figures/PaCMAP/' + hue + '.png',
                    bbox_inches='tight', dpi=300)
    else:
        plt.title("AML Methylome Map, n= " + str(len(score2)),
                  fontsize=12)
        plt.savefig('../Figures/PaCMAP/PaCMAP_Projection.png',
                    bbox_inches='tight', dpi=300)

    return (plt.show())


def plot_confusion_matrix_individual(clf, x_test, y_test, title='Classification results', tick_fontsize=10, label_fontsize=10):

    from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

    sns.set_theme(style='white')
    predictions = clf.predict(x_test)
    cm = confusion_matrix(y_test, predictions, labels=clf.classes_, normalize='true')
    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=clf.classes_)
    disp.plot(cmap='Blues', values_format='.2f', xticks_rotation='vertical', colorbar=False)

    # Decrease the font size of the numbers inside the confusion matrix
    for texts in disp.text_:
        for text in texts:
            text.set_fontsize(label_fontsize)

    # Decrease the font size of the tick labels
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)

    # Increase the size of the plot
    fig = plt.gcf()
    fig.set_size_inches(5, 5)

    # Add title and axis names and place title in the middle
    plt.title(title +', n=' + str(len(x_test)), fontsize=10, fontweight='bold', pad=10, x=-0.2)
    plt.xlabel('Predicted', fontsize=10, fontweight='bold')
    plt.ylabel('Actual', fontsize=10, fontweight='bold')

    # remove x tick labels
    plt.gca().axes.xaxis.set_ticklabels([])

    plt.show()


def plot_confusion_matrix_stacked(clf, x_train, y_train, x_test, y_test, 
                                  title='', 
                                  tick_fontsize=10, label_fontsize=10,
                                  figsize=(10, 5)):

    sns.set_theme(style='white')

    fig, axs = plt.subplots(1, 2, figsize=figsize, sharey=True)

    fig.subplots_adjust(wspace=0.05)  # Adjust the width space

    for i, (ax, x, y, subset) in enumerate(zip(axs, [x_train, x_test], [y_train, y_test],
                                                ['Train with 10-fold CV', 'Validation'])):
        predictions = clf.predict(x)
        print(f'Overall accuracy score in {subset}: {accuracy_score(y, predictions):.3f}')
        cm = confusion_matrix(y, predictions, labels=clf.classes_, normalize='true')
        disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=clf.classes_)
        disp.plot(cmap='Blues', values_format='.2f', xticks_rotation='vertical', colorbar=False, ax=ax)
    

        # Set font size of the numbers inside the confusion matrix
        for texts in disp.text_:
            for text in texts:
                text.set_fontsize(label_fontsize)

        ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
        ax.set_title(subset + ', n=' + str(len(x)), fontsize=label_fontsize+1, pad=10)
        ax.set_xlabel('Predicted label', fontsize=label_fontsize+1)

        # Add y labels
        if i == 0:
            ax.set_ylabel('True label', fontsize=label_fontsize+1)
        else:
            ax.set_ylabel('')

        # remove x tick labels
        ax.xaxis.set_ticklabels([])
        
    plt.show()


def create_color_dict(df, columns, colors):
    """
    Create a color dictionary using unique values from the specified DataFrame columns.
    
    Args:
    df (pd.DataFrame): DataFrame to extract unique values from.
    columns (list): Columns to extract unique values from.
    colors (list): List of colors to assign to unique values.
    
    Returns:
    dict: A dictionary with unique values as keys and corresponding colors as values.
    """
    unique_values = pd.concat([df[col] for col in columns]).unique()
    return {unique_values[i]: colors[i % len(colors)] for i in range(len(unique_values))}

# -*- coding: utf-8 -*-
"""
Produces simple Sankey Diagrams with matplotlib.
@author: Anneya Golob & marcomanz & pierre-sassoulas
Customized by Francisco Marchi
                      .-.
                 .--.(   ).--.
      <-.  .-.-.(.->          )_  .--.
       `-`(     )-'             `)    )
         (o  o  )                `)`-'
        (      )                ,)
        ( ()  )                 )
         `---"\    ,    ,    ,/`
               `--' `--' `--'
                |  |   |   |
                |  |   |   |
                '  |   '   |
"""

from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

class pySankeyException(Exception):
    pass
class NullsInFrame(pySankeyException):
    pass
class LabelMismatch(pySankeyException):
    pass

def check_data_matches_labels(labels, data, side):
    if len(labels >0):
        if isinstance(data, list):
            data = set(data)
        if isinstance(data, pd.Series):
            data = set(data.unique().tolist())
        if isinstance(labels, list):
            labels = set(labels)
        if labels != data:
            msg = "\n"
            if len(labels) <= 20:
                msg = "Labels: " + ",".join(labels) +"\n"
            if len(data) < 20:
                msg += "Data: " + ",".join(data)
            raise LabelMismatch('{0} labels and data do not match.{1}'.format(side, msg))
    

def sankey(left, right, leftWeight=None, rightWeight=None, colorDict=None,
           leftLabels=None, rightLabels=None, aspect=4, rightColor=False,
           fontsize=14, figure_name=None,closePlot=False):
    '''
    Make Sankey Diagram showing flow from left-->right

    Inputs:
        left = NumPy array of object labels on the left of the diagram
        right = NumPy array of corresponding labels on the right of the diagram
            len(right) == len(left)
        leftWeight = NumPy array of weights for each strip starting from the
            left of the diagram, if not specified 1 is assigned
        rightWeight = NumPy array of weights for each strip starting from the
            right of the diagram, if not specified the corresponding leftWeight
            is assigned
        colorDict = Dictionary of colors to use for each label
            {'label':'color'}
        leftLabels = order of the left labels in the diagram
        rightLabels = order of the right labels in the diagram
        aspect = vertical extent of the diagram in units of horizontal extent
        rightColor = If true, each strip in the diagram will be be colored
                    according to its left label
    Ouput:
        None
    '''
    if leftWeight is None:
        leftWeight = []
    if rightWeight is None:
        rightWeight = []
    if leftLabels is None:
        leftLabels = []
    if rightLabels is None:
        rightLabels = []
    # Check weights
    if len(leftWeight) == 0:
        leftWeight = np.ones(len(left))

    if len(rightWeight) == 0:
        rightWeight = leftWeight

    plt.figure()
    plt.rc('text', usetex=False)
    plt.rc('font', family='sans-serif')

    # Create Dataframe
    if isinstance(left, pd.Series):
        left = left.reset_index(drop=True)
    if isinstance(right, pd.Series):
        right = right.reset_index(drop=True)
    df = pd.DataFrame({'left': left, 'right': right, 'leftWeight': leftWeight,
                       'rightWeight': rightWeight}, index=range(len(left)))
    
    if len(df[(df.left.isnull()) | (df.right.isnull())]):
        raise NullsInFrame('Sankey graph does not support null values.')

    # Identify all labels that appear 'left' or 'right'
    allLabels = pd.Series(np.r_[df.left.unique(), df.right.unique()]).unique()

    # Identify left labels
    if len(leftLabels) == 0:
        leftLabels = pd.Series(df.left.unique()).unique()
        leftCounts = df.left.value_counts()  # Count occurrences of each label in 'left'
    else:
        check_data_matches_labels(leftLabels, df['left'], 'left')

    # Identify right labels
    if len(rightLabels) == 0:
        rightLabels = pd.Series(df.right.unique()).unique()
        rightCounts = df.right.value_counts()  # Count occurrences of each label in 'right'
    else:
        check_data_matches_labels(leftLabels, df['right'], 'right')
        
    # If no colorDict given, make one
    if colorDict is None:
        colorDict = {}
        pal = "hls"
        cls = sns.color_palette(pal, len(allLabels))
        for i, l in enumerate(allLabels):
            colorDict[l] = cls[i]
    else:
        missing = [label for label in allLabels if label not in colorDict.keys()]
        if missing:
            raise RuntimeError('colorDict specified but missing values: '
                                '{}'.format(','.join(missing)))

    # Determine widths of individual strips
    ns_l = defaultdict()
    ns_r = defaultdict()
    for l in leftLabels:
        myD_l = {}
        myD_r = {}
        for l2 in rightLabels:
            myD_l[l2] = df[(df.left == l) & (df.right == l2)].leftWeight.sum()
            myD_r[l2] = df[(df.left == l) & (df.right == l2)].rightWeight.sum()
        ns_l[l] = myD_l
        ns_r[l] = myD_r

    # Determine positions of left label patches and total widths
    widths_left = defaultdict()
    for i, l in enumerate(leftLabels):
        myD = {}
        myD['left'] = df[df.left == l].leftWeight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['left']
        else:
            myD['bottom'] = widths_left[leftLabels[i - 1]]['top'] + 0.02 * df.leftWeight.sum()
            myD['top'] = myD['bottom'] + myD['left']
            topEdge = myD['top']
        widths_left[l] = myD

    # Determine positions of right label patches and total widths
    widths_right = defaultdict()
    for i, l in enumerate(rightLabels):
        myD = {}
        myD['right'] = df[df.right == l].rightWeight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['right']
        else:
            myD['bottom'] = widths_right[rightLabels[i - 1]]['top'] + 0.02 * df.rightWeight.sum()
            myD['top'] = myD['bottom'] + myD['right']
            topEdge = myD['top']
        widths_right[l] = myD

    # Total vertical extent of diagram
    xMax = topEdge / aspect

    # Draw vertical bars on left and right of each  label's section & print label
    for l in leftLabels:
        plt.fill_between(
            [-0.02 * xMax, 0],
            2 * [widths_left[l]['bottom']],
            2 * [widths_left[l]['bottom'] + widths_left[l]['left']],
            color=colorDict[l],
            alpha=0.99
        )
        plt.text(
        -0.05 * xMax,
        widths_left[l]['bottom'] + 0.5 * widths_left[l]['left'],
        "{}, n={}".format(l, leftCounts[l]),  # Append sample size to the label
        {'ha': 'right', 'va': 'center'},
        fontsize=fontsize
    )
    for l in rightLabels:
        plt.fill_between(
            [xMax, 1.02 * xMax], 2 * [widths_right[l]['bottom']],
            2 * [widths_right[l]['bottom'] + widths_right[l]['right']],
            color=colorDict[l],
            alpha=0.99
        )
        plt.text(
        1.05 * xMax,
        widths_right[l]['bottom'] + 0.5 * widths_right[l]['right'],
        "{}, n={}".format(l, rightCounts[l]),  # Append sample size to the label
        {'ha': 'left', 'va': 'center'},
        fontsize=fontsize
    )

    # Plot strips
    for l in leftLabels:
        for l2 in rightLabels:
            lc = l
            if rightColor:
                lc = l2
            if len(df[(df.left == l) & (df.right == l2)]) > 0:
                # Create array of y values for each strip, half at left value, half at right, convolve
                ys_d = np.array(50 * [widths_left[l]['bottom']] + 50 * [widths_right[l2]['bottom']])
                ys_d = np.convolve(ys_d, 0.05 * np.ones(20), mode='valid')
                ys_d = np.convolve(ys_d, 0.05 * np.ones(20), mode='valid')
                ys_u = np.array(50 * [widths_left[l]['bottom'] + ns_l[l][l2]] + 50 * [widths_right[l2]['bottom'] + ns_r[l][l2]])
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')

                # Update bottom edges at each label so next strip starts at the right place
                widths_left[l]['bottom'] += ns_l[l][l2]
                widths_right[l2]['bottom'] += ns_r[l][l2]
                plt.fill_between(
                    np.linspace(0, xMax, len(ys_d)), ys_d, ys_u, alpha=0.65,
                    color=colorDict[lc]
                )
    plt.gca().axis('off')
    plt.gcf().set_size_inches(6, 6)
    if figure_name!=None:
        plt.savefig("{}.png".format(figure_name), bbox_inches='tight', dpi=150)
    if closePlot:
        plt.close()


def draw_sankey_plot(df, col1, col2, colors, title, fontsize=10, fig_size=(10,10),
                     column_title=True, column_title_pos = (0.1,0.9)):
    """
    Create a sankey plot using the specified DataFrame columns.
    
    Args:
    df (pd.DataFrame): DataFrame to plot.
    col1, col2 (str): Column names to plot.
    colors (list): List of colors to assign to unique values.
    title (str): Title for the plot.
    fontsize (int): Font size to use for text in the plot. Default is 10.
    fig_size (tuple): Figure size. Default is (10,10).
    column_title (bool): Whether to add column titles to the plot. Default is True.
    column_title_pos (tuple): Relative positions of the column titles. Default is (0.1,0.9).

    """

    # Create color dictionary for unique values in the specified columns
    color_dict = create_color_dict(df, [col1, col2], colors)

    # Fill NaN values with 'Unknown'
    df = df.fillna('Unknown')
    # Add the unknown class
    color_dict['Unknown'] = colors[-1]

    # Create the sankey plot
    sankey(left=df[col1], right=df[col2], aspect=20, colorDict=color_dict, fontsize=fontsize)

    # Get the current figure to add text labels
    fig = plt.gcf()

    if column_title:
        # Add labels at relative positions
        fig.text(column_title_pos[0], 0.88, col1, ha='right', va='center', fontweight='bold', fontsize=fontsize)
        fig.text(column_title_pos[1], 0.88, col2, ha='left', va='center', fontweight='bold', fontsize=fontsize)

    # Set figure size
    fig.set_size_inches(fig_size[0], fig_size[1])

    # Add title
    plt.title(title + ', n=' +
                 str(len(df)), pad=30, fontsize=fontsize)
    
    # Set global font size and family
    plt.rcParams.update({
        'font.size': fontsize,
    })


def get_custom_color_palette():
    list = [
    '#1f77b4',  # Vivid blue
    '#ff7f0e',  # Vivid orange 
    '#2ca02c',  # Vivid green
    '#d62728',  # Vivid red
    '#9467bd',  # Vivid purple 
    '#7f7f7f',  # Medium gray
    '#e377c2',  # Pink
    '#e7ba52',  # Light orange
    '#bcbd22',  # Olive
    '#17becf',  # Light blue
    '#393b79',  # Dark blue
    '#8c564b',  # Brown
    '#f7b6d2',  # Light pink
    '#c49c94',  # Light brown
    '#a2769e',   # Soft purple
    '#dbdb8d',  # Pale yellow
    '#9edae5',  # Pale cyan
    '#c5b0d5',  # Pale purple
    '#c7c7c7',  # Light gray
    '#ff9896',  # Light red
    '#637939',  # Dark olive
    '#aec7e8',  # Light blue
    '#ffbb78',  # Light orange
    '#98df8a',  # Light green
    '#7c231e',  # Dark red
    '#3d6a3d',  # Dark green
    '#f96502',  # Deep orange
    '#6d3f7d',  # Deep purple
    '#6b4423',  # Dark brown
    '#d956a6'   # Hot pink
    ]
    return list