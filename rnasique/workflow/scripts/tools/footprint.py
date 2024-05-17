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
    """
    Compares two reactivity profiles (two 'conditions')

    Args:
        cond1_path: aggreact_tsv file of condition 1
        cond1_name: name of condition 1
        cond2_path: aggreact_tsv file of condition 2
        cond1_name: name of condition 2
        
    Returns:
        footprint (DataFrame): results prepared for plotting by footprint.plot_reactivity()
        footprint_csv (DataFrame): results in table format (to be written to tsv file)
    """

    cols = [
        "seqNum",
        "sequence",
        "pvalue",
        "difference",
        "ratio",
        "significant_higher",
        "significant_lower"
    ]
    indexes = ["seqNum", "sequence"]
    cond1 = pd.read_csv(cond1_path, sep="\t").set_index(["seqNum", "sequence"])
    cond2 = pd.read_csv(cond2_path, sep="\t").set_index(["seqNum", "sequence"])
    
    #joins rows of the two aggreact tables with same seqNum (and nucleotide)
    # (WARNING: if cond2 starts at lower seqNum, these positions are at the end!)
    footprint = pd.concat(
        {cond1_name: cond1, cond2_name: cond2}, names=["condition"], axis=1
    )
    res = list()
    for index, row in footprint.iterrows():
        if row.loc[cond1_name]["desc"] in ("accepted", "warning") and row.loc[
            cond2_name]["desc"] in ("accepted", "warning"):
            # get reactivities of the replicates
            rs1 = list(row.loc[cond1_name].iloc[:row.loc[cond1_name].index.get_loc("mean")])
            rs2 = list(row.loc[cond2_name].iloc[:row.loc[cond2_name].index.get_loc("mean")])
            # compare by ttest only if at least one condition has more than
            # 1 replicates
            if len(rs1)>1 or len(rs2)>1:
                _, pvalue = scipy.stats.ttest_ind(
                    rs1, rs2,
                    equal_var=True,
                    alternative="two-sided",)
            else:
                pvalue = 0 # set to 0 in order to leave decision to other criteria

            difference = row.loc[cond2_name]["mean"] - row.loc[cond1_name]["mean"]
            mean_sum = row.loc[cond2_name]["mean"] + row.loc[cond1_name]["mean"]
            if mean_sum!=0:
                # check again: what is the exact idea of ratio??
                ratio = np.abs(difference) / (mean_sum)
            else: ratio = 0

            significant_higher = 'NO'
            significant_lower = 'NO'
            
            if pvalue < ttest_pvalue_thres and ratio > ratio_thres:
                if difference > diff_thres:
                   significant_higher = 'YES'
                elif difference < -diff_thres:
                   significant_lower = 'YES'

            curres = pd.DataFrame(
                [[index[0], index[1], pvalue, difference, ratio, significant_higher, significant_lower]],
                columns=cols,).set_index(indexes)
            # todo: deprecated; it is inefficient to concat every row
            res.append(curres)
        else:
            curres = pd.DataFrame(
                [[index[0], index[1], np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]],
                columns=cols,).set_index(indexes)
            res.append(curres)

    #res = pd.DataFrame([], columns=cols).set_index(indexes)
    res = pd.concat(res)

    footprint = pd.concat(
        {cond1_name: cond1, cond2_name: cond2, "ttest": res},
        names=["condition"],
        axis=1,)

    footprint_csv = pd.concat(
        {cond1_name: cond1["mean"], cond2_name: cond2["mean"], "analysis": res},
        axis=1)
    footprint_csv = footprint_csv.sort_values(by=["seqNum"])

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
    width = 100, ## sequence positions per row in plot
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

    significant = footprint[[('ttest', 'significant_higher'),('ttest', 'significant_lower')]].reset_index()
    significant = significant.sort_values(by=[('seqNum', '')]).reset_index().drop([('index','')], axis=1)
    significant["xlabel"] = (significant[('seqNum', '')].astype(str) + "\n" + significant[('sequence','')].astype(str))
    unidsignificant = significant.drop([('seqNum', ''), ('sequence','')], axis=1)
   
    difference = footprint[('ttest', 'difference')].reset_index()
    difference = difference.sort_values(by=[('seqNum', '')]).reset_index().drop([('index','')], axis=1)
    difference["xlabel"] = (difference[('seqNum', '')].astype(str) + "\n" + difference[('sequence','')].astype(str))
    difference = difference.drop([('seqNum', ''), ('sequence','')], axis=1)
 
    regions = np.linspace(0, int(len(unidmeans)/width)*width, int(len(unidmeans)/width)+1)
    regions = list(map(int, regions))
    if width/2 > len(unidmeans) - regions[-1] +1:
        regions[-1] = len(unidmeans)
    else:
        regions.append(len(unidmeans))
#    subplot_width = (len(regions)-2) * [1] + [(regions[-1]-regions[-2]+1)/width]


    if len(regions) > 2:
        fig, axes = plt.subplots(len(regions)-1, 1, figsize=(len(footprint) / 3, 4*(len(regions)-1)))
        for j in range(len(regions)-1) :
            axes[j].set_xticks(unidmeans.index[regions[j]:regions[j+1]])
            axes[j].set_xticklabels(unidmeans["xlabel"][regions[j]:regions[j+1]])

            for i in range(len(unidmeans)):
                if i >= regions[j] and i < regions[j+1]:
                    if unidsignificant.loc[i, ('ttest', 'significant_higher')] == 'YES':
                        axes[j].axvspan(unidmeans.index[i] - 0.45, unidmeans.index[i] + 0.45, color='yellow',alpha=0.5)
                    elif unidsignificant.loc[i, ('ttest', 'significant_lower')] == 'YES':
                        axes[j].axvspan(unidmeans.index[i] - 0.45, unidmeans.index[i] + 0.45, color='orange',alpha=0.3)

                    bar_color1 = 'skyblue' if unidmeans.iloc[i, 0] >= 0 else 'lightgrey'
                    bar_color2 = 'royalblue' if unidmeans.iloc[i, 1] >= 0 else 'lightgrey'

                    axes[j].bar(unidmeans.index[i] - 0.225, unidmeans.iloc[i, 0], yerr=unidstdev.iloc[i, 0], capsize=2, width=0.45, color=bar_color1, align='center', label=f"{cond1_name}")
                    axes[j].bar(unidmeans.index[i] + 0.225, unidmeans.iloc[i, 1], yerr=unidstdev.iloc[i, 1], capsize=2, width=0.45, color=bar_color2, align='center', label=f"{cond2_name}")

            axes[j].axhline(y=0.4, color="orange", linestyle="-", label="Medium Reactivity")
            axes[j].axhline(y=0.7, color="red", linestyle="-", label="High reactivity")
            axes[j].axhline(y=0.0, color="silver", linestyle="-")

            legend_elements = [
                Patch(facecolor="skyblue", edgecolor='grey',label=f"{cond1_name}"),
                Patch(facecolor="royalblue", edgecolor='grey',label=f"{cond2_name}"),
                Patch(facecolor="lightgrey", edgecolor='grey',label="Undetermined"),
                Patch(facecolor="yellow", edgecolor='grey',label=f"{cond2_name} sign. higher"),
                Patch(facecolor="orange", edgecolor='grey',label=f"{cond2_name} sign. lower"),
                Line2D([0], [0], color="red", label="High reactivity threshold"),
                Line2D([0], [0], color="orange", label="Medium reactivity threshold"),
            ]
            axes[j].legend(handles=legend_elements, loc="upper left", ncol=14)
            axes[j].set_title(f"Nucleotides {unidmeans['xlabel'][regions[j]].split()[0]} - {unidmeans['xlabel'][regions[j+1]-1].split()[0]}",loc='left')
            axes[j].set_xlim([unidmeans.index[regions[j]] - 1, unidmeans.index[regions[j+1]-1] + 1])
            axes[j].set_ylim([min(min(np.nanmin(unidmeans.iloc[:, 0]-unidstdev.iloc[:, 0]),np.nanmin(unidmeans.iloc[:, 1]-unidstdev.iloc[:, 1]))*1.1,-0.15), \
                max(np.nanmax(unidmeans.iloc[:, 0]+unidstdev.iloc[:, 0]),np.nanmax(unidmeans.iloc[:, 1]+unidstdev.iloc[:, 1]))*1.2])
            ax2 = axes[j].twinx()
            ax2.set_ylim([min(min(np.nanmin(unidmeans.iloc[:, 0]-unidstdev.iloc[:, 0]),np.nanmin(unidmeans.iloc[:, 1]-unidstdev.iloc[:, 1]))*1.1,-0.15), \
                max(np.nanmax(unidmeans.iloc[:, 0]+unidstdev.iloc[:, 0]),np.nanmax(unidmeans.iloc[:, 1]+unidstdev.iloc[:, 1]))*1.2])
            plt.tight_layout()
        plt.suptitle(title)
        plt.savefig(output, format=format)

        fig, axes = plt.subplots(len(regions)-1, 1, figsize=(len(footprint) / 3, 4*(len(regions)-1)))

        for j in range(len(regions)-1) :
            axes[j].set_xticks(unidmeans.index[regions[j]:regions[j+1]])
            axes[j].set_xticklabels(unidmeans["xlabel"][regions[j]:regions[j+1]])

            for i in range(len(unidmeans)):
                if i >= regions[j] and i < regions[j+1]:
                    if unidsignificant.loc[i, ('ttest', 'significant_higher')] == 'YES' or unidsignificant.loc[i, ('ttest', 'significant_lower')] == 'YES':
                        axes[j].bar(unidmeans.index[i], difference.iloc[i, 0], width=0.8, color = 'gold',edgecolor='orange', align='center', label="Significant difference")
                    elif unidsignificant.loc[i, ('ttest', 'significant_higher')] == 'NO' and unidsignificant.loc[i, ('ttest', 'significant_lower')] == 'NO':
                        axes[j].bar(unidmeans.index[i], difference.iloc[i, 0], width=0.8, color = 'palegreen', edgecolor='forestgreen', align='center', label="Difference")
                    else:
                        axes[j].bar(unidmeans.index[i], -0.1, width=0.8, color = 'lightgrey', align='center', label="Undetermined")

            axes[j].axhline(y=0.0, color="silver", linestyle="-")

            legend_elements = [
                Patch(facecolor="palegreen", label="Difference"),
                Patch(facecolor="gold", label="Significant difference"),
                Patch(facecolor="lightgrey", label="Undetermined"),
            ]
            axes[j].legend(handles=legend_elements, loc="upper left")
            axes[j].set_title(f"Nucleotides {unidmeans['xlabel'][regions[j]].split()[0]} - {unidmeans['xlabel'][regions[j+1]-1].split()[0]}",loc='left')
            axes[j].set_xlim([difference.index[regions[j]] - 1, difference.index[regions[j+1]-1] + 1])
            axes[j].set_ylim([min(np.nanmin(difference.iloc[:, 0])*1.1, -0.15), np.nanmax(difference.iloc[:, 0])*1.1])
            ax2 = axes[j].twinx()
            ax2.set_ylim([min(np.nanmin(difference.iloc[:, 0])*1.1, -0.15), np.nanmax(difference.iloc[:, 0])*1.1])
            plt.tight_layout()
    else: # <= 2 regions
        fig, ax = plt.subplots(figsize=(len(footprint) / 3, 4*(len(regions)-1)))
        ax.set_xticks(unidmeans.index[regions[0]:regions[1]])
        ax.set_xticklabels(unidmeans["xlabel"][regions[0]:regions[1]])

        for i in range(len(unidmeans)):
            if unidsignificant.loc[i, ('ttest', 'significant_higher')] == 'YES':
                ax.axvspan(unidmeans.index[i] - 0.45, unidmeans.index[i] + 0.45, color='yellow',alpha=0.5)
            elif unidsignificant.loc[i, ('ttest', 'significant_lower')] == 'YES':
                ax.axvspan(unidmeans.index[i] - 0.45, unidmeans.index[i] + 0.45, color='orange',alpha=0.3)

            bar_color1 = 'skyblue' if unidmeans.iloc[i, 0] >= 0 else 'lightgrey'
            bar_color2 = 'royalblue' if unidmeans.iloc[i, 1] >= 0 else 'lightgrey'

            ax.bar(unidmeans.index[i] - 0.225, unidmeans.iloc[i, 0], yerr=unidstdev.iloc[i, 0], capsize=2, width=0.45, color=bar_color1, align='center', label=f"{cond1_name}")
            ax.bar(unidmeans.index[i] + 0.225, unidmeans.iloc[i, 1], yerr=unidstdev.iloc[i, 1], capsize=2, width=0.45, color=bar_color2, align='center', label=f"{cond2_name}")

            ax.axhline(y=0.4, color="orange", linestyle="-", label="Medium Reactivity")
            ax.axhline(y=0.7, color="red", linestyle="-", label="High reactivity")
            ax.axhline(y=0.0, color="silver", linestyle="-")

        legend_elements = [
            Patch(facecolor="skyblue", edgecolor='grey',label=f"{cond1_name}"),
            Patch(facecolor="royalblue", edgecolor='grey', label=f"{cond2_name}"),
            Patch(facecolor="lightgrey", edgecolor='grey',label="Undetermined"),
            Patch(facecolor="yellow", edgecolor='grey',label=f"{cond2_name} sign. higher"),
            Patch(facecolor="orange", edgecolor='grey',label=f"{cond2_name} sign. lower"),
            Line2D([0], [0], color="red", label="High reactivity threshold"),
            Line2D([0], [0], color="orange", label="Medium reactivity threshold"),
        ]
        ax.legend(handles=legend_elements, loc="upper left",ncol=14)
        ax.set_title(f"Nucleotides {unidmeans['xlabel'][regions[0]].split()[0]} - {unidmeans['xlabel'][regions[1]-1].split()[0]}",loc='left')
        ax.set_xlim([unidmeans.index[regions[0]] - 1, unidmeans.index[regions[1]-1] + 1])
        ax.set_ylim([min(min(np.nanmin(unidmeans.iloc[:, 0]-unidstdev.iloc[:, 0]),np.nanmin(unidmeans.iloc[:, 1]-unidstdev.iloc[:, 1]))*1.1,-0.15), \
            max(np.nanmax(unidmeans.iloc[:, 0]+unidstdev.iloc[:, 0]),np.nanmax(unidmeans.iloc[:, 1]+unidstdev.iloc[:, 1]))*1.2])
        ax2 = ax.twinx()
        ax2.set_ylim([min(min(np.nanmin(unidmeans.iloc[:, 0]-unidstdev.iloc[:, 0]),np.nanmin(unidmeans.iloc[:, 1]-unidstdev.iloc[:, 1]))*1.1,-0.15), \
            max(np.nanmax(unidmeans.iloc[:, 0]+unidstdev.iloc[:, 0]),np.nanmax(unidmeans.iloc[:, 1]+unidstdev.iloc[:, 1]))*1.2])
        plt.tight_layout()
        plt.suptitle(title)
        plt.savefig(output, format=format)

        fig, ax = plt.subplots(figsize=(len(footprint) / 3, 4*(len(regions)-1)))

        ax.set_xticks(unidmeans.index[regions[0]:regions[1]])
        ax.set_xticklabels(unidmeans["xlabel"][regions[0]:regions[1]])

        for i in range(len(unidmeans)):
            if unidsignificant.loc[i, ('ttest', 'significant_higher')] == 'YES' or unidsignificant.loc[i, ('ttest', 'significant_lower')] == 'YES':
                ax.bar(unidmeans.index[i], difference.iloc[i, 0], width=0.8, color = 'gold',edgecolor='orange', align='center', label="Significant difference")
            elif unidsignificant.loc[i, ('ttest', 'significant_higher')] == 'NO' and unidsignificant.loc[i, ('ttest', 'significant_lower')] == 'NO':
                ax.bar(unidmeans.index[i], difference.iloc[i, 0], width=0.8, color = 'palegreen', edgecolor='forestgreen',align='center', label="Difference")
            else:
                ax.bar(unidmeans.index[i], -0.1, width=0.8, color = 'lightgrey', align='center', label="Undetermined")


        ax.axhline(y=0.0, color="silver", linestyle="-")

        legend_elements = [
            Patch(facecolor="palegreen", label="Difference"),
            Patch(facecolor="gold", label="Significant difference"),
            Patch(facecolor="lightgrey", label="Undetermined"),
        ]
        ax.legend(handles=legend_elements, loc="upper left")
        ax.set_title(f"Nucleotides {unidmeans['xlabel'][regions[0]].split()[0]} - {unidmeans['xlabel'][regions[1]-1].split()[0]}",loc='left')
        ax.set_xlim([difference.index[regions[0]] - 1, difference.index[regions[1]-1] + 1])
        ax.set_ylim([min(np.nanmin(difference.iloc[:, 0])*1.1, -0.15), np.nanmax(difference.iloc[:, 0])*1.1])
        ax2 = ax.twinx()
        ax2.set_ylim([min(np.nanmin(difference.iloc[:, 0])*1.1, -0.15), np.nanmax(difference.iloc[:, 0])*1.1])
        plt.tight_layout()
    plt.suptitle(dif_title)
    plt.savefig(diff_output, format=format)


def footprint_2D_plot(infile, higher, lower, outfile, show=True):
    vecs=[higher,lower]
    for i,vec in enumerate(vecs):
        pos=[i for i,x in vec.items() if x=='YES']
        vecs[i] = str(pos).replace(' ','')[1:-1]

    import os 
    conda_prefix = os.environ.get("CONDA_PREFIX")

    varna_cmd = ['java', '-jar', f'{conda_prefix}/lib/varna/VARNA.jar',
                 #'-sequenceDBN', sequence,
                 #'-structureDBN', structure,
                 '-i', infile,
                 '-o', outfile,
                 '-basesStyle1', 'fill=#ff00aa',
                 '-applyBasesStyle1on', vecs[0],
                 '-basesStyle2', 'fill=#aaff00',
                 '-applyBasesStyle2on', vecs[1]]

    print(varna_cmd)
    import subprocess
    subprocess.run(varna_cmd)

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
    structure = None, ## dbn structure file
    structure_plot=None,
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

    if structure_plot is not None:
        print(f"Write structure plot {structure_plot} based on file {structure}.")
        assert(structure is not None)
        footprint = footprint_csv.reset_index().set_index('seqNum')
        higher = footprint['analysis']['significant_higher']
        lower = footprint['analysis']['significant_lower']
        
        footprint_2D_plot(structure, higher, lower, outfile=structure_plot)

def main():
    return fire.Fire(footprint_main)


if __name__ == "__main__":
    main()
