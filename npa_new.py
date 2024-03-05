#!/usr/bin/env python
# coding: utf-8

__author__ = "Arjit M; amisra2@illinois.edu"
__version__ = "Feb 2 2024"

import numpy as np
import pandas as pd
import re
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats
import os
from npa_consts import P29_COMPONENTS, FUP_LOC, FUP_RES, TUMOR_VARS, PREFIX_TO_LABELS, PREFIX_IS_ONE_HOT, \
    PHYS_HLTH_SUMMARY, MENT_HLTH_SUMMARY
from npa_helpers import *

# http://www.healthmeasures.net/media/kunena/attachments/257/PROMIS29_Scoring_08082018.pdf
PAIN_INT_MEAN = 2.31
PAIN_INT_STD = 2.34

df = pd.read_csv("npadata_race.csv")
col_names = df.columns.to_numpy()
fup_cols = [c for c in col_names if re.match("fup_symptoms___\d*", c)]
bsl_cols = [c for c in col_names if re.match("bsl_symptoms___\d*", c)]
fup_total = np.sum(df[fup_cols], axis=1)
fup_total[~np.any(df[fup_cols], axis=1)] = np.nan

bsl_total = np.sum(df[bsl_cols], axis=1)
bsl_total[~np.any(df[bsl_cols], axis=1)] = np.nan

df["fup_total"] = fup_total
df["bsl_total"] = bsl_total

# Convert t_scores to z_scores. Here, T mean is 50 and std is 10.
df['p29_pf_z_score'] = ( df['p29_pf_t_score'] - 50 ) / 10
df['p29_anxiety_z_score'] = ( df['p29_anxiety_t_score'] - 50 ) / 10
df['p29_depression_z_score'] = ( df['p29_depression_t_score'] - 50 ) / 10
df['p29_fatigue_z_score'] = ( df['p29_fatigue_t_score'] - 50 ) / 10
df['p29_sd_z_score'] = ( df['p29_sd_t_score'] - 50 ) / 10
df['p29_social_z_score'] = ( df['p29_social_t_score'] - 50 ) / 10
df['p29_pain_z_score'] = ( df['p29_pain_t_score'] - 50 ) / 10

df['p29_pain_int_z_score'] = (df['p29_global07'] - PAIN_INT_MEAN) / PAIN_INT_STD
df['p29_pain_int_t_score'] = df['p29_pain_int_z_score'] * 10 + 50

df['pain_avg_z'] = np.mean(df[['p29_pain_int_z_score', 'p29_pain_z_score']], axis=1)

df['emotional_dist_z'] = np.mean(df[['p29_depression_z_score', 'p29_anxiety_z_score']], axis=1)

summ_in = df[[
    'p29_pf_z_score',
    'pain_avg_z',
    'p29_social_z_score',
    'p29_fatigue_z_score',
    'p29_sd_z_score',
    'emotional_dist_z',
]]

df["p29_Mental_Health_Summ"] = (MENT_HLTH_SUMMARY @ summ_in.T) * 10 + 50
df["p29_Physical_Health_Summ"] = (PHYS_HLTH_SUMMARY @ summ_in.T) * 10 + 50

col_names = df.columns.to_numpy()

adf = df.groupby("pt_study_id").first()
adf = adf.copy()
adf["symptom_diff"] = adf["bsl_total"] - adf["fup_total"]

adf.to_csv("npa_expanded.csv")

def do_anova(var_head, outputs, out_max, one_hot=True, min_size=15, p_thresh=0.05, result_dir='results_point',
             out_min=None, plot_mode=True):

    num_out = len(outputs)
    outputs = list(outputs)

    if not out_min:
        out_min = [0] * num_out

    categories = []

    if one_hot:
        # 1-hot coded variables have the form category___x with corresponding 0/1 T/F value.
        # e.g. tumor_loc___3 refers to a specific tumor location.
        # Here, isolate categories (__x) for a variable (tumor_loc) for which sufficient data exists, at least
        # MIN_SIZE entries with non-nan outputs
        for col in col_names:
            M = re.match("{0}_*\d+".format(var_head), col)
            if M:
                if len(adf.loc[adf[M.string] == 1][outputs].dropna()) > min_size:
                    categories.append(M.string)

    else:
        # Non 1-hot coded variables have integer values corresponding to categories, such that category = x.
        # e.g. metastatic_no = 5 refers to a specific non-metastatic tumor.
        # Here, isolate categories (x) corresponding to variable (metastatic_no) for which sufficient data exists,
        # at least MIN_SIZE entries with non-nan outputs.
        categories = [C for C, n in zip(*np.unique(adf[var_head], return_counts=True)) if n > min_size]
        categories = [C for C in categories if len(adf.loc[adf[var_head] == C][outputs].dropna()) > min_size]


    cat_num = len(categories)

    if cat_num < 2:
        print("{0} ANOVA not performed. Insuffucient categories.".format(var_head))
        return

    # Physical meaning of integer value.
    labeller = PREFIX_TO_LABELS.get(var_head)

    # If results folder does not exist, create folder in CWD.
    os.makedirs("{0}/{1}".format(result_dir, var_head), exist_ok=True)

    outtxt = []
    outcsv = []
    if plot_mode:
        fig, ax = plt.subplots( (num_out + 1) // 2 , 2, figsize=(19.2, 7 * ((num_out + 1) // 2)), sharex="all")
    plt.yticks()

    lgnd_txt = ""
    outtxt.append(var_head + "\n")
    outtxt.append("=========== Key ===========\n")
    catnum_to_labnum = dict()
    cat_labs = []

    for k, C in enumerate(categories):
        if one_hot:
            lab_num = int(re.match("{0}_+(\d+)".format(var_head), C).group(1))
        else:
            # If dictionary is empty, simply use value in database. Useful for variables whose numerical values are
            # not group numbers e.g. number of lesions.
            lab_num = int(C)

        lab = labeller.get(lab_num, lab_num)
        catnum_to_labnum[k] = lab_num
        cat_labs.append(str(lab))
        outtxt.append("Group {0}: \t  {1}\n".format(k, lab))
        lgnd_txt = lgnd_txt + "Group {0}: {1}\n".format(k, lab)

    outtxt.append('\n')
    labnum_to_catnum = {v:k for k, v in catnum_to_labnum.items()}

    # Initialize first row of output csv with blank column, then names of categories as columns
    outcsv.append("\t".join([''] + cat_labs + ["ANOVA"]))

    '''
    Map for formatting. Given 3-4 plots, want label to be at height = 1/6 of figure.
    1, 2 -> 2;
    3, 4 -> 3;
    5, 6 -> 4;
    '''
    if plot_mode:
        bot_algn  = lambda k: (k + 1) // 2 + 1
        plt.subplots_adjust(bottom = 1 / bot_algn(num_out))
        fig.text(0.5, 0.5 / bot_algn(num_out), lgnd_txt, transform=fig.transFigure, ha="center", va='center', ma='left',
                 bbox=dict(ec="black", fill=False))

    for out_i, output in enumerate(outputs):

        outcsv_row = [output]
        to = []
        for C in categories:
            if one_hot:
                to.append(adf.loc[adf[C] == 1][output].dropna().to_numpy())
            else:
                to.append(adf.loc[adf[var_head] == C][output].dropna().to_numpy())

        # ANOVA with all groups.
        f, p = stats.f_oneway(*to)
        outtxt.append("\n\n##################################################\n")
        outtxt.append(output + "\n\n")
        outtxt.append("p = {0}\nf = {1}\n".format(p, f))
        outtxt.append("\n")

        diff_str = ['a' for _ in range(cat_num)]
        if p < p_thresh:
            res = stats.tukey_hsd(*to)
            Pij, T = res.pvalue, res.statistic
            outtxt.append("=========== P Values Tukey HSD ===========\n")
            outtxt.append("_________________ All ____________________\n")
            for i in range(cat_num):
                for j in range(i + 1, cat_num):
                    outtxt.append("({0}, {1}): {2}\n".format(i, j, str(Pij[i][j])))
            outtxt.append('\n')

            outtxt.append("_____________ Significant ________________\n")
            for i in range(cat_num):
                for j in range(i + 1, cat_num):
                    if Pij[i][j] < p_thresh:
                        outtxt.append("({0}, {1})\np: {2}\nt: {3}\n".format(i, j, str(Pij[i][j]), str(T[i][j])))

            outtxt.append('\n')

            G = Pij < p_thresh
            diff_str = getGroupLabels(G)

            for k, C in enumerate(categories):
                if one_hot:
                    lab = labeller.get(int(re.match("{0}_+(\d+)".format(var_head), C).group(1)))
                else:
                    lab = labeller.get(int(C), int(C))
                outtxt.append("{0}^({1})\nN = {2}\n".format(lab, diff_str[k], len(to[k])))
            outtxt.append("\n")

        else:
            pass

        outtxt.append('=========== Summary ===========\n')
        for k, v in enumerate(to):
            outtxt.append("Group: {0}\nMean: {1}\nStd: {2}\nN: {3}\n".format(k, np.mean(v), np.std(v), len(v)))
            outcsv_row.append("{0} ({1}) N={2}".format(round(np.mean(v), 2), diff_str[k], len(v)))

        outtxt.append('\n')
        if p < p_thresh:
            outcsv_row.append("p={0} f={1}".format(round(p, 4), round(f, 2)))
        else:
            outcsv_row.append("p>{0}".format(p_thresh))
        outcsv.append("\t".join(outcsv_row))

        if plot_mode:
            make_ind_plots(num_out, ax, out_i, out_max, out_min, output, to, catnum_to_labnum, labnum_to_catnum,
                           labeller, result_dir, var_head, diff_str)
        # End for each output

    if plot_mode:
        fig.suptitle(var_head)
        fig.savefig("{0}/{1}.jpg".format(result_dir, var_head), dpi=600)
        plt.close()

    with open("{0}/{1}_anova_tHSD.txt".format(result_dir, var_head), 'w') as outfile:
        for ln in outtxt:
            outfile.write(ln)

    with open("{0}/{1}_anova_tHSD.tsv".format(result_dir, var_head), 'w') as outfile:
        for ln in outcsv:
            outfile.write(ln + "\n")


for var in PREFIX_TO_LABELS.keys():
    print(var)
    do_anova(
        var,
        P29_COMPONENTS.get("OUTPUTS"),
        P29_COMPONENTS.get("OUTPUT_MAX"),
        result_dir="results_raw_outputs",
        one_hot=PREFIX_IS_ONE_HOT.get(var)
    )

    do_anova(
        var,
        P29_COMPONENTS.get("OUTPUTS_t"),
        [100] * len(P29_COMPONENTS.get("OUTPUTS_t")),
        result_dir="results_t_outputs",
        one_hot=PREFIX_IS_ONE_HOT.get(var)
    )

    do_anova(
        var,
        ["p29_Mental_Health_Summ", "p29_Physical_Health_Summ"],
        [100,100],
        result_dir="results_summary",
        one_hot=PREFIX_IS_ONE_HOT.get(var)
    )

    do_anova(
        var,
        ["fup_total", "bsl_total", "symptom_diff"],
        out_max=[10, 10, 5],
        one_hot=PREFIX_IS_ONE_HOT.get(var),
        result_dir="results_symptoms",
        out_min=[0, 0, -5],
    )





