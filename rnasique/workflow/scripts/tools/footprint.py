#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
import scipy.stats
import fire

# def ratio_sig_test(footprint, ttest_pvalue_thres=0.05, diff_thres=0.2, ratio_thres=0.2):
#     footprint = pd.DataFrame(footprint)
#     footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
#         lambda row: row["delta"]
#         if row["pvalue"] <= ttest_pvalue_thres
#         and np.abs(row["delta"]) >= diff_thres
#         and row["ratio"] >= ratio_thres
#         else np.NaN,
#         axis=1,
#     )
#     return footprint


# def zfactor_sig_test(footprint, ttest_pvalue_thres=0.01, zfactor_thres=0):
#     footprint = pd.DataFrame(footprint)
#     footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
#         lambda row: row["delta"]
#         if row["pvalue"] <= ttest_pvalue_thres and row["zfactor"] > zfactor_thres
#         else np.NaN,
#         axis=1,
#     )
#     return footprint


# def ttest_only_sig_test(footprint, ttest_pvalue_thres=0.01):
#     footprint = pd.DataFrame(footprint)
#     footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
#         lambda row: row["delta"] if row["pvalue"] <= ttest_pvalue_thres else np.NaN,
#         axis=1,
#     )
#     return footprint


# def ratio_zfactor_sig_test(
#     footprint,
#     ttest_pvalue_thres=0.01,
#     diff_thres=0.2,
#     ratio_thres=0.2,
#     zfactor_thres=0,
# ):
#     footprint = pd.DataFrame(footprint)
#     footprint.loc[:, ("ttest", "significant_delta")] = footprint["ttest"].apply(
#         lambda row: row["delta"]
#         if row["pvalue"] <= ttest_pvalue_thres
#         and row["zfactor"] > zfactor_thres
#         and np.abs(row["delta"]) >= diff_thres
#         and row["ratio"] >= ratio_thres
#         else np.NaN,
#         axis=1,
#     )
#     return footprint


def footprint_ttest(
    cond1_path,
    cond1_name,
    cond2_path,
    cond2_name,
    deviation_type="stdev",
    # zfactor_nsigma=3,
    ttest_pvalue_thres=0.05,
    diff_thres=0.2,
    ratio_thres=0.2,
):
    cols = [
        "seqNum",
        "sequence",
        "pvalue",
        "difference",
        "ratio",
        "significant",
    ]
    indexes = ["seqNum", "sequence"]
    cond1 = pd.read_csv(cond1_path, sep="\t").set_index(["seqNum", "sequence"])
    cond2 = pd.read_csv(cond2_path, sep="\t").set_index(["seqNum", "sequence"])
    footprint = pd.concat(
        {cond1_name: cond1, cond2_name: cond2}, names=["condition"], axis=1
    )
    res = pd.DataFrame([], columns=cols).set_index(indexes)
    for index, row in footprint.iterrows():
        if row.loc[cond1_name]["desc"] in ("accepted", "warning") and row.loc[
            cond2_name]["desc"] in ("accepted", "warning"):
            stat, pvalue = scipy.stats.ttest_ind(
                list(row.loc[cond1_name].iloc[:row.loc[cond1_name].index.get_loc("mean")]),
                list(row.loc[cond2_name].iloc[:row.loc[cond1_name].index.get_loc("mean")]),
                equal_var=True,
                alternative="two-sided",)
            difference = np.abs(row.loc[cond2_name]["mean"] - row.loc[cond1_name]["mean"])
            ratio = difference / (row.loc[cond2_name]["mean"] + row.loc[cond1_name]["mean"])

            if pvalue < ttest_pvalue_thres and difference >= diff_thres and ratio >= ratio_thres:
                significant = 'YES'
            else:
                significant = 'NO'
            curres = pd.DataFrame(
                [[index[0], index[1], pvalue, difference, ratio, significant]],
                columns=cols,).set_index(indexes)
            res = pd.concat([res, curres])
        else:
            curres = pd.DataFrame(
                [[index[0], index[1], np.NaN, np.NaN, np.NaN, np.NaN]],
                columns=cols,).set_index(indexes)
            res = pd.concat([res, curres])
    footprint = pd.concat(
        {cond1_name: cond1, cond2_name: cond2, "ttest": res},
        names=["condition"],
        axis=1,)

    footprint_csv = pd.concat(
        {cond1_name: cond1["mean"], cond2_name: cond2["mean"], "analysis": res},
        names=["condition"],
        axis=1,)
    footprint_csv = footprint_csv.sort_values(by=["seqNum"]).reset_index()

    return footprint, footprint_csv


def plot_reactivity(
    footprint: pd.DataFrame,
    title="Footprint",
    dif_title="Footprint",
    output="fig.svg",
    diff_output="diff.svg",
    format="svg",
    cond1_name="Condition1",
    cond2_name="Condition2",
):

    replicates = footprint.drop("ttest", axis=1)
    means = replicates.xs("mean", level=1, axis=1).reset_index()
    means = means.sort_values(by=["seqNum"]).reset_index().drop(["index"], axis=1)
    means["xlabel"] = (means["seqNum"].astype(str) + "\n" + means["sequence"].astype(str))
    means = means.replace(-10, -0.1)
    unidmeans = means.drop(["seqNum", "sequence"], axis=1)

    stdev = replicates.xs("stdev", level=1, axis=1).reset_index()
    stdev = stdev.sort_values(by=["seqNum"]).reset_index().drop(["index"], axis=1)
    stdev["xlabel"] = (stdev["seqNum"].astype(str) + "\n" + stdev["sequence"].astype(str))
    unidstdev = stdev.drop(["seqNum", "sequence"], axis=1)

    significant = footprint[('ttest', 'significant')].reset_index()
    significant = significant.sort_values(by=[('seqNum', '')]).reset_index().drop([('index','')], axis=1)
    significant["xlabel"] = (significant[('seqNum', '')].astype(str) + "\n" + significant[('sequence','')].astype(str))
    unidsignificant = significant.drop([('seqNum', ''), ('sequence','')], axis=1)
   
    difference = footprint[('ttest', 'difference')].reset_index()
    difference = difference.sort_values(by=[('seqNum', '')]).reset_index().drop([('index','')], axis=1)
    difference["xlabel"] = (difference[('seqNum', '')].astype(str) + "\n" + difference[('sequence','')].astype(str))
    difference = difference.drop([('seqNum', ''), ('sequence','')], axis=1)
 
    regions = np.linspace(0, int(len(unidmeans)/100)*100, int(len(unidmeans)/100)+1)
    regions = list(map(int, regions))
    if 50 > len(unidmeans) - regions[-1] +1:
        regions[-1] = len(unidmeans)
    else:
        regions.append(len(unidmeans))
#    subplot_width = (len(regions)-2) * [1] + [(regions[-1]-regions[-2]+1)/100]

    fig, axes = plt.subplots(len(regions)-1, 1, figsize=(len(footprint) / 3, 4*(len(regions)-1)))

    for j in range(len(regions)-1) :
        axes[j].set_xticks(unidmeans.index[regions[j]:regions[j+1]])
        axes[j].set_xticklabels(unidmeans["xlabel"][regions[j]:regions[j+1]])
    
        for i in range(len(unidmeans)):
            if i >= regions[j] and i < regions[j+1]:
                if unidsignificant.loc[i, ('ttest', 'significant')] == 'YES':
                    axes[j].axvspan(unidmeans.index[i] - 0.45, unidmeans.index[i] + 0.45, color='gold',alpha=0.3)

                bar_color1 = 'skyblue' if unidmeans.iloc[i, 0] >= 0 else 'lightgrey'
                bar_color2 = 'royalblue' if unidmeans.iloc[i, 1] >= 0 else 'lightgrey'
                
                axes[j].bar(unidmeans.index[i] - 0.225, unidmeans.iloc[i, 0], yerr=unidstdev.iloc[i, 0], capsize=2, width=0.45, color=bar_color1, align='center', label=f"{cond1_name}")
                axes[j].bar(unidmeans.index[i] + 0.225, unidmeans.iloc[i, 1], yerr=unidstdev.iloc[i, 1], capsize=2, width=0.45, color=bar_color2, align='center', label=f"{cond2_name}")

        axes[j].axhline(y=0.4, color="orange", linestyle="-", label="Medium Reactivity")
        axes[j].axhline(y=0.7, color="red", linestyle="-", label="High reactivity")
        axes[j].axhline(y=0.0, color="silver", linestyle="-")

        legend_elements = [
            Patch(facecolor="skyblue", label=f"{cond1_name}", alpha=0.5),
            Patch(facecolor="royalblue", label=f"{cond2_name}", alpha=0.5),
            Patch(facecolor="lightgrey", label="Undetermined", alpha=0.5),
            Patch(facecolor="gold", label="Significant difference",alpha=0.5),
            Line2D([0], [0], color="red", label="High reactivity threshold"),
            Line2D([0], [0], color="orange", label="Medium reactivity threshold"),
        ]
        axes[j].legend(handles=legend_elements, loc="upper left")
        axes[j].set_title(f"Nucleotides {unidmeans['xlabel'][regions[j]].split()[0]} - {unidmeans['xlabel'][regions[j+1]-1].split()[0]}",loc='left')
        axes[j].set_xlim([unidmeans.index[regions[j]] - 1, unidmeans.index[regions[j+1]-1] + 1])
        axes[j].set_ylim([min(min(np.nanmin(unidmeans.iloc[:, 0]-unidstdev.iloc[:, 0]),np.nanmin(unidmeans.iloc[:, 0]-unidstdev.iloc[:, 1]))*1.1,0), \
            max(np.nanmax(unidmeans.iloc[:, 0]+unidstdev.iloc[:, 0]),np.nanmax(unidmeans.iloc[:, 0]+unidstdev.iloc[:, 1]))*1.1])
        ax2 = axes[j].twinx()
        ax2.set_ylim([min(min(np.nanmin(unidmeans.iloc[:, 0]-unidstdev.iloc[:, 0]),np.nanmin(unidmeans.iloc[:, 0]-unidstdev.iloc[:, 1]))*1.1,0), \
            max(np.nanmax(unidmeans.iloc[:, 0]+unidstdev.iloc[:, 0]),np.nanmax(unidmeans.iloc[:, 0]+unidstdev.iloc[:, 1]))*1.1])
        plt.tight_layout()
    plt.suptitle(title)
    plt.savefig(output, format=format)

    fig, axes = plt.subplots(len(regions)-1, 1, figsize=(len(footprint) / 3, 4*(len(regions)-1)))

    for j in range(len(regions)-1) :
        axes[j].set_xticks(unidmeans.index[regions[j]:regions[j+1]])
        axes[j].set_xticklabels(unidmeans["xlabel"][regions[j]:regions[j+1]])

        for i in range(len(unidmeans)):
            if i >= regions[j] and i < regions[j+1]:
                if unidsignificant.loc[i, ('ttest', 'significant')] == 'YES':
                    axes[j].bar(unidmeans.index[i], difference.iloc[i, 0], width=0.5, color = 'blue', align='center', label="Significant difference")
                elif unidsignificant.loc[i, ('ttest', 'significant')] == 'NO':
                    axes[j].bar(unidmeans.index[i], difference.iloc[i, 0], width=0.5, color = 'skyblue', align='center', label="Difference")
                else:
                    axes[j].bar(unidmeans.index[i], -0.1, width=0.5, color = 'lightgrey', align='center', label="Undetermined")

        axes[j].axhline(y=0.0, color="silver", linestyle="-")

        legend_elements = [
            Patch(facecolor="skyblue", label="Difference"),
            Patch(facecolor="blue", label="Significant difference"),
            Patch(facecolor="lightgrey", label="Undetermined"),
        ]
        axes[j].legend(handles=legend_elements, loc="upper left")
        axes[j].set_title(f"Nucleotides {unidmeans['xlabel'][regions[j]].split()[0]} - {unidmeans['xlabel'][regions[j+1]-1].split()[0]}",loc='left')
        axes[j].set_xlim([difference.index[regions[j]] - 1, difference.index[regions[j+1]-1] + 1])
        axes[j].set_ylim([-0.15, np.nanmax(difference.iloc[:, 0])*1.1])
        ax2 = axes[j].twinx()
        ax2.set_ylim([-0.15, np.nanmax(difference.iloc[:, 0])*1.1])
        plt.tight_layout()
    plt.suptitle(dif_title)
    plt.savefig(diff_output, format=format)



def footprint_main(
    cond1: str,
    cond2: str,
    cond1_name: str = None,
    cond2_name: str = None,
    # deviation_type: str = "stdev",
    ttest_pvalue_thres=0.05,
    diff_thres=0.2,
    ratio_thres=0.2,
    output: str = None,
    plot: str = None,
    diff_plot: str = None,
    diff_plot_title="Footprint",
    plot_title="Footprint",
    plot_format="svg",
):
    cond1_name = cond1_name if cond1_name is not None else cond1
    cond2_name = cond2_name if cond2_name is not None else cond2
    # assert deviation_type in ["stdev", "sem", "mad"]

    footprint, footprint_csv = footprint_ttest(
        cond1,
        cond1_name,
        cond2,
        cond2_name,
        # deviation_type=deviation_type,
        ttest_pvalue_thres=0.05,
        diff_thres=0.2,
        ratio_thres=0.2,
    )

    if output is not None:
        footprint_csv.to_csv(output, sep="\t")

    if plot is not None or diff_plot is not None:
        plot_reactivity(
            footprint,
            plot_title,
            diff_plot_title,
            plot,
            diff_plot,
            plot_format,
            # deviation_type=deviation_type,
            cond1_name=cond1_name,
            cond2_name=cond2_name,
        )



def main():
    return fire.Fire(footprint_main)


if __name__ == "__main__":
    main()
