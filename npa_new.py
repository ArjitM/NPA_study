#!/usr/bin/env python
# coding: utf-8

__author__ = "Arjit M; amisra2@illinois.edu"
__version__ = "Jan 14 2024"

import numpy as np
import pandas as pd
import re
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats
from npa_consts import P29_COMPONENTS, FUP_LOC, FUP_RES, FUP_SYMPTOM_LABELS, TUMOR_VARS, TUMOR_LOC_LABELS, \
    APPROACH_LABELS, PREOP_TX_LABELS, INPATIENT_LABELS, DISHARGE_DISP_LABELS, REASON_READMIT_LABELS, PREFIX_TO_LABELS

df = pd.read_csv("npadata_race.csv")
col_names = df.columns.to_numpy()
fup_cols = [c for c in col_names if re.match("fup_symptoms___\d*", c)]

fup_total = np.sum(df[fup_cols], axis=1)
fup_total[~np.any(df[fup_cols], axis=1)] = np.nan

df["fup_total"] = fup_total

MAX_VAL = 20
qol_out_cols = []
for output in P29_COMPONENTS.get("OUTPUTS"):
    qol_out_cols.extend(P29_COMPONENTS.get(output))
    qol_out_cols.append(output)
# comp_df = df[qol_out_cols].dropna(how='any')

p29_total = np.sum(df[P29_COMPONENTS.get("OUTPUTS_POS")], axis=1) \
            + MAX_VAL * len(P29_COMPONENTS.get("OUTPUTS_NEG")) \
            - np.sum(df[P29_COMPONENTS.get("OUTPUTS_NEG")], axis=1)

# All components of survey must be filled out for inclusion
p29_total[np.any(np.isnan(df[qol_out_cols]), axis=1)] = np.nan

df["p29_total"] = p29_total
col_names = df.columns.to_numpy()

tdf = df[TUMOR_VARS + ["pt_study_id"]].groupby("pt_study_id").first()

tsize = tdf[["tsize_diam1","tsize_diam2", "tsize_diam3"]].convert_dtypes(convert_string=True)
tunit = tdf[["tsize_axial_unit_2", "tsize_axial_unit_4", "tsize_axial_unit_3"]]
tunit.dropna(how="all")

tsize.dropna(how='all')

tdf[[
    "tsize_axial", 
    "tsize_coronal", 
    "tsize_sagittal",
#     "tsize_diam1",
#     "tsize_diam2", 
#     "tsize_diam3",
]].dropna(how='all')

# tvol = []  # 3D tumor volume calculated as d1.d2.d3
# dmax = []  # maximum diameter along any axis
# pat3 = "\d+.*\d*\s*x\s*\d+.*\d*s*x\s*\d+.*\d*"  # regex to match a.b x c.d x e.f; a x c.d x e allowed, space agnostic
# pat2 = "\d+.*\d*\s*x\s*\d+.*\d*"                # regex to match a.b x c.d
# x = 0
# y = 0
# z = 0
# for row in tsize.itertuples(name=None, index=False):
#     try:
#         frow = [float(x) for x in row]
#         tvol.append(np.prod(frow))
#         dmax.append(np.max(frow))
#         z = z+1
#     except:
#         x = x + 1
#         d3 = np.sum([bool(re.match(pat3, str(c))) for c in row])
#         #d2 = np.sum([bool(re.match(pat2, str(c))) for c in row])
        
#         if d3 >= 1:
#             y = y+1
#             V = np.max([ np.prod([float(n) for n in re.split("[A-z|\s]+", str(x)) if n]) for x in row])
#             L = []
#             [L.extend([float(n) for n in re.split("[A-z|\s]+", str(x)) if n]) for x in row]
#             M = np.max(L)
#             tvol.append(V)
#             dmax.append(M)
#             print("**", row, V, M)
#         else:
#             tvol.append(np.nan)
#             dmax.append(np.nan)
#             print(row)

tvol = []  # 3D tumor volume calculated as d1.d2.d3
dmax = []  # maximum diameter along any axis

for row in tdf.itertuples(index=False):
    try:
        frow = [float(d) for d in [row.tsize_diam1, row.tsize_diam2, row.tsize_diam3]]
        # convert to mm, assume default is mm for nan values. key: 1 = cm; 2 = mm
        units_pow = np.sum([u==1 for u in [row.tsize_axial_unit_2, row.tsize_axial_unit_4, row.tsize_axial_unit_3]])
        tvol.append(np.prod(frow) * 10**units_pow)
        dmax.append(np.max(frow))
    except:
        # print([row.tsize_diam1, row.tsize_diam2, row.tsize_diam3])
        tvol.append(np.nan)
        dmax.append(np.nan)


tdf["tvol"] = tvol
tdf["dmax"] = dmax

adf = df.groupby("pt_study_id").first()

# COL_ORDER = {v: k for k,v in enumerate(col_names)}
# dt = df[sorted(set(col_names) - set(tumor_vars), key=COL_ORDER.get)].merge(tdf, on="pt_study_id", how='inner')

def do_anova(var_head, min_size=30):

    categories = []
    for col in col_names:
        M = re.match("{0}_*\d+".format(var_head), col)
        if M:
            if len(adf.loc[adf[M.string] == 1]) > min_size:
                categories.append(M.string)
    cat_num = len(categories)

    with open("results/{0}_anova_tHSD.txt".format(var_head), 'w') as outfile:

        for output in P29_COMPONENTS.get("OUTPUTS"):
            to = []

            for C in categories:
                to.append(adf.loc[adf[C] == 1][output].dropna().to_numpy())

            _, p = stats.f_oneway(*to)
            outfile.write("\n\n##################################################\n")
            outfile.write(output + "\np = " + str(p))
            outfile.write("\n")

            if p < 1e-3:

                Pij = stats.tukey_hsd(*to).pvalue
                print("=========== P Values Tukey HSD ===========")
                for i in range(cat_num):
                    for j in range(i + 1, cat_num):
                        outfile.write("({0}, {1}): {2}\n".format(i+1, j+1, str(Pij[i][j])))
                outfile.write('\n')

            outfile.write('=========== Summary ===========\n')
            [outfile.write("Group: {0}\nMean: {1}\nStd: {2}\n".format(i, np.mean(v), np.std(v))) for i, v in enumerate(to)]
            outfile.write('\n')

            label_data = {}
            labeller = PREFIX_TO_LABELS.get(var_head)
            for C, data in zip(categories, to):
                num = int(re.match("{0}_+(\d+)".format(var_head), C).group(1))
                label_data[labeller.get(num)] = data

            plt.figure()
            sns.boxenplot(data=label_data, color='black')
            plt.savefig("results/{0}_{1}.jpg".format(var_head, output))


do_anova("tumor_laterality")




