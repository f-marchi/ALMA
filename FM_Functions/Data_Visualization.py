"""
This module implements functions for several commonly used data visualization techniques.

"""

# Import libraries

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def draw_kaplan_meier(scorename, df, save_plot=False,
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
    f, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 10))

    # Define survival curve categories
    groups = df[scorename + ' Categorical']
    ix = (groups == 'High')

    # Fit the Kaplan Meier to each category of groups (kmf1 and kmf2)
    def surv_curves(i, t, e):

        T = df[t]
        E = df[e]
        kmf1 = KaplanMeierFitter()
        kmf1.fit(T[~ix], E[~ix], label='Low-' + scorename + ', n=' +
                 str(len(df[df[scorename + ' Categorical'].isin(['Low'])])))
        ax = kmf1.plot_survival_function(
            ax=i, show_censors=True, ci_show=show_ci)

        kmf2 = KaplanMeierFitter()
        kmf2.fit(T[ix], E[ix], label='High-' + scorename + ', n=' +
                 str(len(df[df[scorename + ' Categorical'].isin(['High'])])))
        ax = kmf2.plot_survival_function(
            ax=i, show_censors=True, ci_show=show_ci)

        # Calculate Hazard Ratio (HZ) and p-value (p)
        X_CPH = df[[scorename + '_cat_bin', t, e]]
        cph = CoxPHFitter()
        HZ = cph.fit(X_CPH, t, event_col=e)
        hz = HZ.hazard_ratios_[0]
        p = HZ.summary['p'][0]

        # Annotate HZ and p
        i.annotate(f'Hazard Ratio: {hz:.4f}\np-value: {p:.4f}',
                   xy=(9.75, 0.085), xycoords='data',
                   ha='right', va='center', fontsize=11,
                   bbox={'boxstyle': 'round', 'facecolor': 'none',
                         'edgecolor': 'lightgray'})

        # Add risk counts below the graph
        if add_risk_counts == True:
            add_at_risk_counts(kmf1, kmf2, ax=i)

        # Save Survival Function Table
        if save_survival_table == True:
            surv1 = kmf1.survival_function_.join(kmf1.confidence_interval_)
            surv2 = kmf2.survival_function_.join(kmf2.confidence_interval_)
            surv3 = surv1.join(surv2, how='outer')
            surv3.to_csv('../Figures/Kaplan_Meiers/KM_OS_SurvivalTable_' +
                         scorename + '_' + trialname + '_' + str(len(df)) + '.csv')

        i.set_ylim(0, 1)
        i.set_ylabel("est. probability of survival $\hat{S}(t)$")

    surv_curves(i=ax[0], t='efs.time', e='efs.evnt')
    surv_curves(i=ax[1], t='os.time', e='os.evnt')

    ax[0].set_title('Event-Free Survival', loc='left',
                    pad=10, fontweight='bold')
    ax[1].set_title('Overall Survival', loc='left', pad=10, fontweight='bold')
    # Define Plot Specs
    plt.subplots_adjust(wspace=0, hspace=0.2)
    plt.suptitle("Kaplan-Meiers of " + scorename + " in " + trialname + ", n=" + str(len(df)),
                 fontsize='medium', y=0.94,
                 fontweight='bold')
    plt.xlim(0, 10)
    plt.xlabel("time $t$ (years)")

    # Save plot figure
    if save_plot == True:
        plt.savefig('../Figures/Kaplan_Meiers/' + scorename + '_' + trialname + '_' + str(len(df)) + '.png',
                    bbox_inches='tight', dpi=300)

    return (plt.show())


def draw_forest_plot(time, event, df, save_plot=False, trialname=None, scorename=None):
    """
    Generates a custom forest plot.

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
    scorename: str
        Name of your model.

    Returns:
    --------
        A magnificent forest plot.

    """
    import myforestplot as mfp
    from tableone import TableOne
    import statsmodels.formula.api as smf
    import numpy as np

    fp = df[[scorename+' Categorical',
             'MRD 1 Status',
             'Risk Group',
             'FLT3 ITD',
             'Leucocyte counts (10⁹/L)',
             'Age group (years)',
             time, event]].rename(columns={scorename + ' Categorical': scorename})

    event2 = event.replace('.', '_')
    time2 = time.replace('.', '_')

    if event[0] == 'o':
        event3 = 'Overall-Survival'
        event4 = 'OS'
    else:
        event3 = 'Event-Free Survival'
        event4 = 'EFS'

    fp2 = fp.rename(columns={event: event2,
                             time: time2,
                             'MRD 1 Status': 'MRD_1_Status',
                             'FLT3 ITD': 'FLT3_ITD',
                             'Risk Group': 'Risk_Group',
                             'Leucocyte counts (10⁹/L)': 'WBC_count',
                             'Age group (years)': 'Age_group'})

    res = smf.phreg(formula=time2 + " ~ C("+scorename+",Treatment(reference='Low')) + MRD_1_Status + C(Risk_Group,Treatment(reference='Low Risk')) + FLT3_ITD + WBC_count + Age_group",
                    data=fp2, status=event2).fit()

    res2 = res.summary(xname=[scorename+'-High',
                              'MRD 1 Status-Positive',
                              'Risk Group-High Risk',
                              'Risk Group-Standard Risk',
                              'FLT3 ITD-Yes',
                              'Leucocyte counts (10⁹/L)-≥30',
                              'Age group (years)-≥10']).tables[1]

    res3 = res2.set_index(res2.index.str.split(pat='-', expand=True))

    mytable = TableOne(data=fp.drop(columns=[event, time]),
                       pval=False, missing=True, overall=True,
                       label_suffix=False, order={scorename: ['High'],
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
        plt.savefig('../Figures/Forest_Plots/' + scorename + '_' + trialname + '_' + str(len(df)) + '_' + event4 + '.png',
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
        ax.set_xlabel('Diagnostic samples, n='+str(df.shape[1]), fontsize=10)
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
