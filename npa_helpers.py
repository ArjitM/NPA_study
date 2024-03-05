import numpy as np
import pandas as pd
from npa_consts import TUMOR_VARS
import seaborn as sns
import matplotlib as mpl
from matplotlib import pyplot as plt
from textwrap import wrap

cseq = mpl.color_sequences.get('tab20b')[0::5] + mpl.color_sequences.get('tab20c')[0::5] + \
       mpl.color_sequences.get('tab20b')[2::5] + mpl.color_sequences.get('tab20c')[2::5] + \
       mpl.color_sequences.get('tab20b')[1::5] + mpl.color_sequences.get('tab20c')[1::5] + \
       mpl.color_sequences.get('tab20b')[3::5] + mpl.color_sequences.get('tab20c')[3::5]

def getGroupLabels(G):
    eq_groups = set()
    for i in range(len(G)):
        # for each category i, make an equivalent group with categories that do not significantly differ
        # from category i.
        eq_groups.add(tuple(np.where(~G[i])[0]))

    for i in range(len(G)):
        for j in np.where(G[i])[0]:

            del_groups = set()
            add_grps = set()

            for group in eq_groups.union(add_grps):

                # If categories i and j are significantly different, must split all equivalent groups where
                # both i and j are present.
                # for e.g., consider a group {A B C D}. If p(B, C) < thresh, then  split the group into
                # new groups {A B D} and {A C D}.
                if i in group and j in group:
                    del_groups.add(group)

                    # First new group after split. Add iff non-empty
                    ng = tuple(x for x in group if x != i)
                    if ng:
                        add_grps.add(ng)

                    # Second new group after split. Add iff non-empty
                    ng = tuple(x for x in group if x != j)
                    if ng:
                        add_grps.add(ng)

            for grp in del_groups:
                eq_groups.remove(grp)

            for ag in add_grps:
                eq_groups.add(ag)

    # After adding new groups, some groups may be redundant. If so, retain only the greatest superset.
    # e.g., if equivalent groups include {1, 3, 4} and {3, 4} then {3, 4} is redundant because 3 and 4
    # are already labelled as equivalent.
    redundant = set()
    for grp1 in eq_groups:
        for grp2 in eq_groups:
            if set(grp1) < set(grp2):
                redundant.add(grp1)

    eq_groups = eq_groups - redundant

    cat_str = [list() for _ in range(len(G))]
    for k, group in enumerate(sorted(eq_groups)):
        for i in group:
            cat_str[i].append(chr(ord('a') + k))

    diff_str = [''.join(sorted(x)) for x in cat_str]
    return diff_str


def generateTumorDf(df):

    tdf = df[TUMOR_VARS + ["pt_study_id"]].groupby("pt_study_id").first()
    tsize = tdf[["tsize_diam1","tsize_diam2", "tsize_diam3"]].convert_dtypes(convert_string=True)
    tunit = tdf[["tsize_axial_unit_2", "tsize_axial_unit_4", "tsize_axial_unit_3"]]
    tunit.dropna(how="all")
    tsize.dropna(how='all')

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
            tvol.append(np.nan)
            dmax.append(np.nan)


    tdf["tvol"] = tvol
    tdf["dmax"] = dmax

    tdf.to_csv("npa_tumor_data.csv")


def make_ind_plots(num_out, ax, out_i, out_max, out_min, output, to, catnum_to_labnum, labnum_to_catnum, labeller,
                   result_dir, var_head, diff_str):
    if num_out > 2:
        ax_plt = ax[out_i // 2][out_i % 2]
    elif num_out == 2:
        ax_plt = ax[out_i]
    else:
        ax_plt = ax

    cpal = ['silver' for _ in diff_str]
    unq_lbls = np.unique(diff_str)
    if len(unq_lbls) > 1:
        lmap = {x: y for x, y in zip(sorted(unq_lbls), cseq)}
        cpal = [lmap.get(g) for g in diff_str]

    sns.pointplot(ax=ax_plt, data={k: v for k, v in enumerate(to)},
                  color='black', marker='D', markersize=1.2,
                  # palette=cpal,
                  estimator='mean', errorbar=("ci", 95), linestyle='none', capsize=0.6,
                  err_kws={'linewidth': 0.5, 'solid_capstyle': 'butt', 'color': 'black'})

    bp = sns.barplot(ax=ax_plt, data={k: v for k, v in enumerate(to)}, errorbar=None, alpha=0)

    # for i, patch in enumerate(bp.patches):
    #     bp.text(
    #         patch.get_x() + patch.get_width() / 2,
    #         patch.get_height() + out_max[out_i] / 10,
    #         "\n".join(wrap(diff_str[i], 3)),
    #         ha='center', va='bottom',
    #         )

    ax_plt.set_title(output)
    ax_plt.set_ylim(ymin=out_min[out_i], ymax=out_max[out_i])

    fig_ind = plt.figure(figsize=(10, 12))
    ax_ind = fig_ind.add_subplot()

    """
    Individual plots ordered by size.
    """

    # 0: Category_i values; 1: Category_i+1 values; ...
    data = {k: v for k, v in enumerate(to)}

    # A: 5, B: 3, C: 6, ... such that mean(A) > mean(B) > ...
    size_map = {chr(ord('A')+k): catnum_to_labnum.get(v) for k, v in enumerate(
        sorted(data.keys(), key=lambda m: np.mean(data.get(m)), reverse=True)
    )}

    # A: Category_i+5 values, ...
    size_order_data = {k: data.get(labnum_to_catnum.get(size_map.get(k))) for k in size_map.keys()}

    grp_lbls = [diff_str[labnum_to_catnum.get(lab_i)] for lab_i in size_map.values()]
    # cpal = mpl.color_sequences.get('tab20')[:len(np.unique(grp_lbls))]
    # for i, gl in enumerate()

    cpal = ['silver' for _ in grp_lbls]
    unq_lbls = np.unique(grp_lbls)
    if len(unq_lbls) > 1:
        lmap = {x: y for x, y in zip(sorted(unq_lbls), cseq)}
        cpal = [lmap.get(g) for g in grp_lbls]
        # cpal = ['orange' if len(g) == 1 else 'lightblue' for g in grp_lbls]


    fig_ind.subplots_adjust(bottom=0.5)

    # sns.pointplot(ax=ax_ind, data=size_order_data,
    #               color='black', marker='d', markersize=1,
    #               estimator='mean', errorbar=("ci", 95), linestyle='none', capsize=0.6,
    #               err_kws={'linewidth': 0.6, 'solid_capstyle': 'butt', 'color': 'gray', 'alpha': 0.5})

    B = sns.barplot(ax=ax_ind, data=size_order_data, alpha=1, estimator='mean',
                    palette=cpal,
                    errorbar=("ci", 95), capsize=0.6,
                    err_kws={'linewidth': 0.6, 'solid_capstyle': 'butt', 'color': 'black', 'alpha': 0.8},
                    )
    ax_ind.set_title("{0} by {1}".format(output, var_head))
    ax_ind.set_ylim(ymin=out_min[out_i], ymax=out_max[out_i])


    lgnd_ind = ""
    for k, v in size_map.items():
        lgnd_ind = lgnd_ind + "Group {0}: {1}\n".format(k, labeller.get(v, v))

    fig_ind.text(0.5, 0.05, lgnd_ind, transform=fig_ind.transFigure, ha="center", va='bottom',
             ma='left',
             bbox=dict(ec="black", fill=False))

    # for i, p in enumerate(B.containers):
    #     print(i)
    for i, patch in enumerate(B.patches):
        # plt.bar_label(B.containers[i], fontsize=5, fmt='%.1f', padding=20)
        B.text(
            patch.get_x() + patch.get_width() / 2,
            patch.get_height() + out_max[out_i] / 10,
            "\n".join(wrap(grp_lbls[i], 3)),
            ha='center', va='bottom',
        )

    fig_ind.savefig("{0}/{1}/{2}.jpg".format(result_dir, var_head, output), dpi=600)

