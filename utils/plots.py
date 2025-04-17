"""
This module implements custom plotting functions for data visualization.
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
from bokeh.plotting import figure, show, output_notebook, gridplot
from bokeh.models import Legend, LegendItem
import scipy.stats as stats
from scipy.stats import fisher_exact, chi2_contingency
from sklearn.metrics import (confusion_matrix, ConfusionMatrixDisplay, accuracy_score,
                             precision_recall_fscore_support, roc_auc_score, cohen_kappa_score,
                               roc_curve, auc)


def draw_kaplan_meier(df,
                      model_name,
                      trialname=None,
                      plot_efs=True,
                      plot_os=True,
                      figsize=(8, 10),
                      add_risk_counts=False,
                      save_survival_table=False,
                      show_ci=False,
                      suptitle=0.94,
                      status1= 'High',
                      status2= 'Low',):
    """
    Plot Kaplan-Meier curves for EFS, OS, or both. Preserves all original functionality.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the survival data.
    model_name : str
        Name of the column in df used to define 'High'/'Low' groups.
    trialname : str, optional
        Name of the trial (for figure titles/filenames).
    plot_efs : bool, optional
        Whether to plot EFS curves.
    plot_os : bool, optional
        Whether to plot OS curves.
    save_plot : bool, optional
        Whether to save the resulting plot.
    figsize : tuple, optional
        Size of the figure.
    add_risk_counts : bool, optional
        Whether to add risk tables.
    save_survival_table : bool, optional
        Whether to save the survival tables as CSV.
    show_ci : bool, optional
        Whether to show confidence intervals in the Kaplan-Meier curves.
    suptitle : float, optional
        Y-coordinate of the main title.
    status1 : str, optional
        Label for the 'High' or 'positive' group in the Kaplan-Meier plot.
    status2 : str, optional
            Label for the 'Low' or "negative" group in the Kaplan-Meier plot.
    """

    sns.set_theme(style='white')
    # Import libraries for Kaplan Meier
    from lifelines.plotting import add_at_risk_counts
    from lifelines import KaplanMeierFitter
    from lifelines import CoxPHFitter

    # Collect survival endpoints to plot
    survival_info = []
    if plot_efs:
        survival_info.append(('efs.time at 5y', 'efs.evnt at 5y', 'Event-Free Survival'))
    if plot_os:
        survival_info.append(('os.time at 5y', 'os.evnt at 5y', 'Overall Survival'))

    # If no endpoints specified, return
    if not survival_info:
        return

    # Set up subplots
    fig, axes = plt.subplots(nrows=len(survival_info), ncols=1,
                             sharex=True, figsize=figsize)
    # If only one subplot, axes is not a list
    if len(survival_info) == 1:
        axes = [axes]

    groups = df[model_name]
    ix_high = (groups == status1)

    # Cox model requires numeric indicator
    df[model_name + '_int'] = df[model_name].map({status1: 1, status2: 0})

    def plot_surv(axis, time_col, event_col, title):
        kmf_low = KaplanMeierFitter()
        kmf_high = KaplanMeierFitter()

        T = df[time_col]
        E = df[event_col]

        kmf_low.fit(T[~ix_high],
                    E[~ix_high],
                    label= status2 + ' ' + model_name + ', n=' + str((df[model_name] == status2).sum()))
        kmf_high.fit(T[ix_high],
                     E[ix_high],
                     label=status1 + ' ' + model_name + ', n=' + str((df[model_name] == status1).sum()))

        kmf_low.plot_survival_function(ax=axis, show_censors=True, ci_show=show_ci)
        kmf_high.plot_survival_function(ax=axis, show_censors=True, ci_show=show_ci)

        # CoxPHFitter for hazard ratio and p-value
        try:
            X_cph = df[[model_name + '_int', time_col, event_col]].dropna()
            cph = CoxPHFitter()
            cph.fit(X_cph, time_col, event_col=event_col)
            hz = cph.hazard_ratios_[model_name + '_int']
            p = cph.summary.loc[model_name + '_int', 'p']
            ci_lower = cph.summary.loc[model_name + '_int', 'exp(coef) lower 95%']
            ci_upper = cph.summary.loc[model_name + '_int', 'exp(coef) upper 95%']
            axis.annotate(
                f'HR: {hz:.4f} (95% CI: {ci_lower:.2f}, {ci_upper:.2f})\n p-value: {p:.4f}',
                xy=(4.90, 0.1),
                xycoords='data',
                ha='right', va='center', fontsize=11,
                bbox={'boxstyle': 'round', 'facecolor': 'none', 'edgecolor': 'lightgray'})
        except:
            pass

        # Add risk counts
        if add_risk_counts:
            add_at_risk_counts(kmf_low, kmf_high, ax=axis)

        # Save survival table
        if save_survival_table:
            surv_low = kmf_low.survival_function_.join(kmf_low.confidence_interval_)
            surv_high = kmf_high.survival_function_.join(kmf_high.confidence_interval_)
            merged = surv_low.join(surv_high, how='outer', rsuffix=f'_{status1}')
            merged.to_csv(f'KM_{time_col}_SurvivalTable_'
                          f'{model_name}_{trialname}_{len(df)}.csv')

        axis.set_title(title, loc='left', pad=10, fontweight='bold', fontsize=10)
        axis.set_ylim(0, 1)
        axis.set_xlim(0, 5)
        axis.set_ylabel("est. probability of survival $\hat{S}(t)$")

    # Plot each requested endpoint
    for (ax_i, (time_col, event_col, ttl)) in zip(axes, survival_info):
        plot_surv(ax_i, time_col, event_col, ttl)

    # Label the bottom axis
    axes[-1].set_xlabel("time $t$ (years)")

    # Main title
    plt.suptitle(model_name + " in " + str(trialname) + ", n=" + str(len(df)),
                 fontsize=11, y=suptitle, fontweight='bold')
    plt.subplots_adjust(wspace=0, hspace=0.2)

    return plt.show()


def draw_forest_plot(time, event, df, trialname=None, model_name=None):
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
             'BM leukemic blasts (%)', 
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
                             'BM leukemic blasts (%)': 'BM_blasts',
                             'Age group (years)': 'Age_group'})

    res = smf.phreg(formula=time2 + " ~ C("+model_name+",Treatment(reference='Low')) + C(MRD_1_Status) + C(Risk_Group,Treatment(reference='Low Risk')) + C(FLT3_ITD) + C(WBC_count) + C(BM_blasts) + C(Age_group)",
                    data=fp2, status=event2).fit()

    res2 = res.summary(xname=[model_name+'-High',
                              'MRD 1 Status-Positive',
                              'Risk Group-High Risk',
                              'Risk Group-Standard Risk',
                              'FLT3 ITD-Yes',
                              'Leucocyte counts (10⁹/L)-≥30',
                              'BM leukemic blasts (%)->50',
                              'Age group (years)-≥10']).tables[1]

    res3 = res2.set_index(res2.index.str.split(pat='-', expand=True))

    mytable = TableOne(data=fp.drop(columns=[event, time]),
                       pval=False, missing=True, overall=True,
                       label_suffix=False, order={model_name: ['High'],
                                                  'MRD 1 Status': ['Positive'],
                                                  'Risk Group': ['High Risk', 'Standard Risk'],
                                                  'FLT3 ITD': ['Yes'],
                                                  'Leucocyte counts (10⁹/L)': ['≥30'],
                                                  'BM leukemic blasts (%)': ['>50'],
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
    fp.axd[2].set_xlim([1, 10.5])
    fp.axd[2].set_xticks([0, 2, 4, 6, 8, 10])
    fp.axd[2].set_xticklabels(labels=[0, 2, 4, 6, 8, 10], fontdict={'fontsize': 8})
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

    return (plt.show())


def draw_boxplot(df, x, y, order, trialname, hue=None, palette=None, figsize=None):

    sns.set_theme(style='white')
    plt.subplots(figsize=figsize)

    if order == 'auto':
        order2 = df[x].value_counts().index.tolist()
    else:
        order2 = order.copy()

    ax = sns.boxplot(y=y, x=x, data=df,
                     whis=[0, 100], width=.6, orient='v', order=order2, color='white')

    sns.stripplot(y=y, x=x, data=df, size=4, linewidth=0, orient='v',
                  order=order2, hue=hue, palette=palette or 'bright', dodge=True)

    ax.xaxis.grid(False)

    sns.despine(trim=True, left=True)
    plt.title(f'{y} by {x} in {trialname}', fontsize='medium', y=1, fontweight='bold')

    if hue:
        plt.legend(loc='center right', bbox_to_anchor=(1.25, 0.5))

    return plt.show()


def draw_stacked_barplot(df, x, y, order, trialname, hue=None, palette=None, figsize=None, fontsize=10, annotate_p=True):

    sns.set_theme(style='white')
    plt.subplots(figsize=figsize)

    if order == 'auto':
        order2 = df[x].value_counts().index.tolist()
    else:
        order2 = order.copy()

    hue_counts = df.groupby([x, hue])[y].count().unstack(fill_value=0)
    hue_props = hue_counts.divide(hue_counts.sum(axis=1), axis=0)

    colors = palette if palette else sns.color_palette('bright', n_colors=hue_props.shape[1])
    hue_props.loc[order2, :].plot(kind='bar', stacked=True, ax=plt.gca(), color=colors)

    for rect in plt.gca().patches:
        height = rect.get_height()
        if height > 0:
            plt.gca().text(rect.get_x() + rect.get_width() / 2,
                           rect.get_y() + height / 2,
                           '{:.1f}%'.format(height * 100),
                           ha='center', va='center', color='white', fontsize=fontsize)

    plt.gca().yaxis.set_major_formatter(mticker.PercentFormatter(1))
    plt.gca().xaxis.grid(False)
    plt.xticks(rotation=0, ha='center')
    plt.title(f'{y} by {x}{trialname}, n={len(df)}', y=1)
    plt.legend(loc='center right', bbox_to_anchor=(1.25, 0.5))

    # Add p-value and test name if binary
    if annotate_p and hue_counts.shape == (2, 2):
        table = hue_counts.loc[order2].values
        if (table < 5).any():
            stat_test = "Fisher's exact"
            _, p = fisher_exact(table)
        else:
            stat_test = 'Chi-squared'
            _, p, _, _ = chi2_contingency(table)

        plt.text(0.4, -0.3,
                 f'{stat_test}, p={p:.4f}',
                 ha='left', va='center', transform=plt.gca().transAxes, fontsize=fontsize)


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


def plot_confusion_matrix_stacked(
    df_dict,
    true_col,
    pred_col,
    class_col,
    title='',
    tick_fontsize=10,
    label_fontsize=10,
    figsize_per_plot=(5, 3),
):
    """
    Plots confusion matrices for 2 or 3 datasets side by side with true labels accompanied by counts.
    """

    def compute_metrics(y_true, y_pred, is_binary):
        metrics = {'Accuracy': accuracy_score(y_true, y_pred)}
        if is_binary:
            precision, recall, f1, _ = precision_recall_fscore_support(y_true, y_pred, average='binary')
            tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
            specificity = tn / (tn + fp)
            auc_roc = roc_auc_score(y_true, y_pred)
            metrics.update({
                'Sensitivity': recall,
                'Specificity': specificity,
                'Precision': precision,
                'F1-score': f1,
                'AUC-ROC': auc_roc
            })
        else:
            weighted_f1 = precision_recall_fscore_support(y_true, y_pred, average='weighted')[2]
            kappa = cohen_kappa_score(y_true, y_pred)
            metrics.update({
                'Weighted F1': weighted_f1,
                "Cohen's Kappa": kappa
            })
        return metrics

    all_labels = set()
    for df in df_dict.values():
        all_labels.update(df[class_col].dropna().unique())
        all_labels.update(df[pred_col].dropna().unique())
    all_classes = sorted([cls for cls in all_labels if pd.notna(cls)])

    is_binary = len(all_classes) == 2
    display_labels = ["Alive", "Dead"] if is_binary else all_classes

    class_to_int = {cls: i for i, cls in enumerate(all_classes)}

    n_plots = len(df_dict)
    fig, axs = plt.subplots(1, n_plots, figsize=(figsize_per_plot[0] * n_plots, figsize_per_plot[1]))
    if n_plots == 1:
        axs = [axs]

    metrics_dict = {}
    for ax, (subset, df) in zip(axs, df_dict.items()):
        y_true, y_pred = df[true_col].map(class_to_int), df[pred_col].map(class_to_int)
        metrics = compute_metrics(y_true, y_pred, is_binary)
        metrics_dict[subset] = metrics

        cm = confusion_matrix(y_true, y_pred, labels=range(len(all_classes)), normalize='true')
        disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=display_labels)
        disp.plot(cmap='Blues', values_format='.2f', xticks_rotation='vertical', 
                  colorbar=False, ax=ax)

        label_counts = df[true_col].map(lambda x: display_labels[x] if isinstance(x, int) else x).value_counts().to_dict()
        y_tick_labels = [f'{label}, n={label_counts.get(label, 0)}' for label in display_labels]
        ax.set_yticklabels(y_tick_labels, fontsize=label_fontsize)

        for texts in disp.text_:
            for text in texts:
                text.set_fontsize(label_fontsize)

        ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
        ax.set_title(f'{subset}, n={len(df)}', fontsize=label_fontsize+1, pad=10)
        ax.set_xlabel('Predicted label', fontsize=label_fontsize+1)
        ax.set_ylabel('True label' if ax == axs[0] else '', fontsize=label_fontsize+1)
        ax.set_xticklabels([])

    if title:
        fig.suptitle(title, fontsize=label_fontsize+2)

    plt.tight_layout()
    plt.show()

    metrics_df = pd.DataFrame(metrics_dict).T
    metrics_df = metrics_df.applymap(lambda x: f"{x:.3f}" if isinstance(x, (int, float)) else x)
    print("\nMetrics:")
    print(metrics_df.to_markdown())


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

        # -*- coding: utf-8 -*-
    """
    Produces simple Sankey Diagrams with matplotlib.
    @author: Anneya Golob & marcomanz & pierre-sassoulas

    ###Customized by Francisco Marchi###
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
                     column_title=True, column_title_pos = (0.1,0.9), nan_action='drop'):
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
    nan_action (str): What to do with NaN values. Options are 'drop', 'include', or 'keep only'. Default is 'drop'.

    """
    df = df.copy()
    
    # Create color dictionary for unique values in the specified columns
    color_dict = create_color_dict(df, [col1, col2], colors)

    if nan_action == 'drop1':
        df = df.dropna(subset=[col1])
        df[col2] = df[col2].fillna('Not confident')

    elif nan_action == 'drop2':
        df = df.dropna(subset=[col1,col2])

    elif nan_action == 'include':
        df[col1] = df[col1].fillna('Unknown')
        if col2 == 'ALMA Subtype':
            df[col2] = df[col2].fillna('Not confident')
        else:
            df[col2] = df[col2].fillna('Unknown')

    elif nan_action == 'keep-only':
        df = df[df[col1].isna()]
        # Fill NaN values with 'Unknown'
        df[col1] = df[col1].fillna('Unknown')
        if col2 == 'ALMA Subtype':
            df[col2] = df[col2].fillna('Not confident')
        else:
            df[col2] = df[col2].fillna('Unknown')
    else:
        raise ValueError('Invalid value for nan_action. Options are "drop1", "drop2" ,"include", or "keep only".')

    # Add the unknown class
    color_dict['Unknown'] = colors[-1]
    color_dict['Not confident'] = colors[-1]

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


def plot_roc_auc(df, target, title=None, color_option='colors1'):
    """
    Plots ROC AUC flexibly using Bokeh, with legends sorted by ascending AUC.
    """
    if color_option == 'colors1':
        colors = ['red', 'green', 'blue', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'black']
    elif color_option == 'colors2':
        colors = ['green', 'blue', 'red']
    else:
        colors = ['green', 'red', 'blue']
    
    if title:
        title_ = title + ', n=' + str(len(df))
    else:
        title_ = ''

    p = figure(title=title_,
               x_axis_label='False Positive Rate',
               y_axis_label='True Positive Rate',
               width=425, height=425,
               tools='save,reset,pan')

    p.line([0, 1], [0, 1], line_dash="dashed", color="gray", line_width=1)

    # Calculate AUCs for each column and store them with their associated column name
    aucs = []
    for column in df.columns.difference([target]):
        fpr, tpr, _ = roc_curve(df[target], df[column])
        roc_auc = auc(fpr, tpr)
        aucs.append((roc_auc, column, fpr, tpr))

    # Sort by AUC in ascending order
    aucs.sort(reverse=True)

    # Plot each ROC curve using sorted AUCs
    for (roc_auc, column, fpr, tpr), color in zip(aucs, colors):
        p.line(fpr, tpr, legend_label=f"{column} ({roc_auc:.2f})",
               color=color, line_width=2, alpha=0.8)

    p.legend.location = "bottom_right"
    p.legend.click_policy = "hide"
    p.toolbar.logo = None
    p.legend.label_text_font_size = '8pt'
    p.legend.spacing = 2
    p.xaxis.axis_label_text_font_style = "normal"
    p.yaxis.axis_label_text_font_style = "normal"
    p.legend.background_fill_alpha = 0.8
    p.title.text_font_size = '10pt'

    return p


def get_custom_color_palette():
    return [
        '#ff7f0e', '#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#7f7f7f',
        '#e377c2', '#e7ba52', '#bcbd22', '#17becf', '#393b79', '#8c564b',
        '#f7b6d2', '#c49c94', '#a2769e', '#dbdb8d', '#9edae5', '#c5b0d5',
        '#c7c7c7', '#ff9896', '#637939', '#aec7e8', '#ffbb78', '#98df8a',
        '#7c231e', '#3d6a3d', '#f96502', '#6d3f7d', '#6b4423', '#d956a6'
    ]


def plot_multiclass_roc_auc(df, target_columns, title=None, width=350, height=800, fontsize='10pt'):
    """
    Plots ROC AUC for multiple classes using Bokeh, targeting multiple classes,
    with legend below the plot.
    """

    colors = get_custom_color_palette()

    if title:
        title_ = title + ', n=' + str(len(df))
    else:
        title_ = ''

    # Initialize figure without a legend inside it
    p = figure(title=title_,
               x_axis_label='False Positive Rate',
               y_axis_label='True Positive Rate',
               width=width, height=height,
               tools='save,reset,pan')

    p.line([0, 1], [0, 1], line_dash="dashed", color="gray", line_width=1)

    legend_items = []  # To hold the legend items

    for target, color in zip(target_columns, colors):
        column = 'P(' + target + ')'
        if column in df.columns:  # Check if the probability column exists
            fpr, tpr, _ = roc_curve(df[target], df[column])
            roc_auc = auc(fpr, tpr)
            line = p.line(fpr, tpr,
                          color=color, line_width=2, alpha=0.8)
            # Create legend items manually
            legend_items.append(LegendItem(label=f"{column} ({roc_auc:.2f})", renderers=[line]))

    # Create a Legend with the collected items and add it below the plot
    legend = Legend(items=legend_items, spacing=0, location="left",
                    label_text_font_size=fontsize, click_policy="hide",
                    background_fill_alpha=0.8)
    p.add_layout(legend, 'below')


    p.toolbar.logo = None
    p.xaxis.axis_label_text_font_style = "normal"
    p.yaxis.axis_label_text_font_style = "normal"
    p.xaxis.axis_label_text_font_size = fontsize
    p.yaxis.axis_label_text_font_size = fontsize
    p.title.text_font_size = fontsize

    return p


def draw_scatter_pearson(df, x, y, s, x_label, y_label):

    sns.set_theme(style="white")
    f, ax = plt.subplots(figsize=(3, 3))
    sns.despine(f, left=False, bottom=False)

    # Set font size for all elements
    sns.set_context("notebook", font_scale=0.8)

    # Define scatterplot

    sns.scatterplot(x=x, y=y,
                    palette='flare', s=s,
                    sizes=(1, 8), linewidth=0, alpha=0.8,
                    data=df, ax=ax)

    # Calculate Pearson’s correlation coefficient (r)
    # and its two-tailed p-value (p)

    r,p = stats.pearsonr(df[x],df[y])

    # Replace plt.annotate() with ax.text() for automatic placement within the plot
    ax.text(0.99, 0.03, fr'$\rho$: {r:.4f}'+f'\np-value: {p:.4f}',
        transform=ax.transAxes,
        horizontalalignment='right', verticalalignment='bottom',
        bbox={'boxstyle': 'round', 'facecolor': 'none', 'edgecolor': 'lightgray' })

    # Define plot specs

    plt.xlabel(x + x_label, fontsize = 10)
    plt.ylabel(y + y_label, fontsize = 10)
    plt.title(r"Pearson's correlation ($\rho$) in " + str(len(df)) + " samples",
               fontsize = 10) 

    return(plt.show())