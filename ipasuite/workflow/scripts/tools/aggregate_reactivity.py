#!/usr/bin/env python3
import os
import pandas as pd
import fire
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mp
from matplotlib import colormaps as cm
import itertools
import logging


import seaborn as sns
sns.set_style("dark")

from scipy.stats import linregress

#plt.style.library["seaborn-dark"]

ddof = 1


class ReactivityThreshold:
    INVALID = -0.3

 #   LOW = 0.5
 #   MEDIUM = 0.4
 #   HIGH = 0.7

    COLOR_INVALID = "grey"
    COLOR_NONE = "white"
    COLOR_LOW = "black"
    COLOR_MEDIUM = "yellow"
    COLOR_HIGH = "red"


def prepare_sequence_axis(ax, sequence):
    font_prop = mp.font_manager.FontProperties(
        family="monospace", style="normal", weight="bold", size="4.5"
    )
    for i in range(sequence):
        nuc = sequence[i]
        if nuc == "T":
            nuc = "U"
        color_dict = {"A": "#f20000", "U": "#f28f00", "G": "#00509d", "C": "#00c200"}
        if nuc in color_dict:
            col = color_dict[nuc]
        elif nuc.upper() in color_dict:
            col = color_dict[nuc.upper()]
        else:
            col = "black"
        ax.annotate(
            nuc,
            xy=(i + 1, -0.67),
            fontproperties=font_prop,
            color=col,
            annotation_clip=False,
            horizontalalignment="center",
        )


def plot_aggregation_infos(aggregated: pd.DataFrame, ax):
    first_idx = aggregated.first_valid_index()[0]
    last_idx = aggregated.last_valid_index()[0]
    #    unit = (56 / 4.0) / (last_idx - first_idx)
    fr = 0
    fone = 0
    fnev = 0
    fnc = 0

    for index, row in aggregated.loc[first_idx:last_idx].iterrows():
        idx = index[0]
        if row["desc"] == "warning":
            ax.axvspan(
                (idx - first_idx) - 0.5,
                (idx - first_idx + 0.4),
                alpha=0.3,
                color="orange",
                label="_" * fr + "Warning",
            )
            fr += 1
        if row["desc"] == "one-value-available":
            ax.axvspan(
                (idx - first_idx) - 0.5,
                (idx - first_idx + 0.4),
                alpha=0.3,
                color="yellow",
                label="_" * fone + "Only one value",
            )
            fone += 1
        if row["desc"] == "undetermined":
            ax.axvspan(
                (idx - first_idx) - 0.5,
                (idx - first_idx + 0.4),
                alpha=0.3,
                color="blue",
                label="_" * fnev + "Undetermined",
            )
            fnev += 1
        if row["desc"] == "non-consistant":
            ax.axvspan(
                (idx - first_idx) - 0.5,
                (idx - first_idx + 0.4),
                alpha=0.3,
                color="red",
                label="_" * fnc + "Non consistant values",
            )
            fnc += 1


def plot_aggregate(
    aggregated: pd.DataFrame,
    reactivity_medium,
    reactivity_high,
    fulloutput="fig.full.svg",
    output="fig.svg",
    title: str = "Aggregated reactivity",
    format="svg",
    ):
    aggregated = copy.deepcopy(aggregated)
    aggregated["xlabel"] = (
        aggregated.index.get_level_values("seqNum").astype(str)
        + "\n"
        + aggregated.index.get_level_values("sequence").astype(str)
    )
    replicates = aggregated.loc[
        :,
        aggregated.columns.drop(["used_values", "mean", "stdev", "desc", "sem"]),
    ].replace(-10, np.NaN)

    aggregated = aggregated.sort_index()
    meanstdev = aggregated.loc[:, ["xlabel", "mean", "stdev"]].replace(-10, np.NaN)

    ax = replicates.plot(
        x="xlabel",
        kind="bar",
        width=0.7,
        stacked=False,
        figsize=(len(aggregated) / 3.5, 4),
        align="center",
        xticks=np.arange(0, len(aggregated) + 1, 1),
    )
    ax = meanstdev[["mean", "xlabel"]].plot(
        x="xlabel",
        y="mean",
        drawstyle="steps-mid",
        ax=ax,
        colormap=cm["cubehelix"],
        linewidth=0.5,
    )

    plot_aggregation_infos(aggregated, ax)

    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.margins(0)
    plt.title(title, loc="left", fontsize=20)
    plt.legend(loc="upper left", ncol=14)
    # print(meanstdev.index.get_level_values("seqNum") - 1)
    try:
        ax.errorbar(
            range(0, meanstdev.shape[0]),
            meanstdev["mean"],
            yerr=meanstdev["stdev"],
            fmt="",
            color="k",
            ls="none",
            capsize=2,
            linewidth=0.5,
        )
    except Exception:  # as e:
        logging.error("fail to generate error bar")
        # logging.error(e)

    ax.set_xlabel("Sequence", fontsize=15)
    ax.set_ylabel("Aggregated reactivity", fontsize=15)
    try:
        plt.tight_layout()
        plt.savefig(fulloutput, format=format)
    except ValueError as e:
        logging.error(f"Unable to save fullplot: {e}")
        open(fulloutput, "a").close()

    regions = np.linspace(0, int(len(aggregated)/100)*100, int(len(aggregated)/100)+1)
    regions = list(map(int, regions))
    if 50 > len(aggregated) - regions[-1] +1:
        regions[-1] = len(aggregated)
    else:
        regions.append(len(aggregated))

    aggregated["color"] = ReactivityThreshold.COLOR_NONE
    aggregated.loc[
        (aggregated["mean"] > reactivity_high), "color"
    ] = ReactivityThreshold.COLOR_HIGH
    aggregated.loc[
        (
            (aggregated["mean"] <= reactivity_high)
            & (aggregated["mean"] > reactivity_medium)
        ),
        "color",
    ] = ReactivityThreshold.COLOR_MEDIUM
    aggregated.loc[
        (
            (aggregated["mean"] <= reactivity_medium)
            & (aggregated["mean"] > ReactivityThreshold.INVALID)
        ),
        "color",
    ] = ReactivityThreshold.COLOR_LOW
    aggregated.loc[
        (aggregated["mean"] <= ReactivityThreshold.INVALID), "color"
    ] = ReactivityThreshold.COLOR_INVALID

    aggregated.loc[(aggregated["mean"] == -10), "stdev"] = np.NaN
    aggregated.loc[(aggregated["mean"] == -10), "mean"] = np.NaN
    aggregated["xlabel_rot"] = (
        aggregated.index.get_level_values("seqNum").astype(str)
        + " - "
        + aggregated.index.get_level_values("sequence").astype(str)
    )

   # plt.style.use("default")
    if len(regions) > 2:
        fig, axes = plt.subplots(len(regions)-1, 1, figsize=(45, 4*(len(regions)-1)))
        agg_index = [int(x[0]) for x in aggregated.index]
        for j in range(len(regions)-1) :
            axes[j].set_xticks(agg_index[regions[j]:regions[j+1]])
            axes[j].set_xticklabels(aggregated["xlabel"][regions[j]:regions[j+1]])
            for i in range(len(aggregated)):
                if i >= regions[j] and i < regions[j+1]:
                    if aggregated["color"][aggregated.index[i]] == 'grey':
                        axes[j].bar(agg_index[i],-0.1, yerr = aggregated["stdev"][aggregated.index[i]], capsize=2, width=0.8, alpha = 0.5,edgecolor='black', color=aggregated["color"][aggregated.index[i]], align='center')
                    else:
                        axes[j].bar(agg_index[i], aggregated["mean"][aggregated.index[i]], yerr = aggregated["stdev"][aggregated.index[i]], capsize=2, width=0.8,edgecolor='grey', color=aggregated["color"][aggregated.index[i]], align='center')
            axes[j].axhline(y=0.0, color="silver", linestyle="-")
            legend_elements = [
            Patch(facecolor=ReactivityThreshold.COLOR_HIGH, edgecolor='grey', label="High Reactivity"),
            Patch(facecolor=ReactivityThreshold.COLOR_MEDIUM,edgecolor='grey', label="Medium Reactivity"),
            Patch(facecolor=ReactivityThreshold.COLOR_LOW,edgecolor='grey', label="Low Reactivity"),
            Patch(facecolor=ReactivityThreshold.COLOR_INVALID, alpha = 0.5,edgecolor='black',label="Undetermined"),
            ]
            axes[j].legend(handles=legend_elements, loc="upper left", ncol=16)
            axes[j].set_title(f"Nucleotides {aggregated['xlabel'][aggregated.index[regions[j]]].split()[0]} - {aggregated['xlabel'][aggregated.index[regions[j+1]-1]].split()[0]}",loc='left', fontsize=20)
            axes[j].set_xlabel("Sequence", fontsize=15)
            axes[j].set_ylabel("Average reactivity", fontsize=15)
            axes[j].set_xlim([agg_index[regions[j]] - 1, agg_index[regions[j+1]-1] + 1])
            axes[j].set_ylim([min(np.nanmin(aggregated["mean"]-aggregated["stdev"])*1.1,-0.15), np.nanmax(aggregated["mean"]+aggregated["stdev"])*1.2])
            ax2 = axes[j].twinx()
            ax2.set_ylim([min(np.nanmin(aggregated["mean"]-aggregated["stdev"])*1.1,-0.15), np.nanmax(aggregated["mean"]+aggregated["stdev"])*1.2])
    else:
        fig, ax = plt.subplots(figsize=(45, 4*(len(regions)-1)))
        agg_index = [int(x[0]) for x in aggregated.index]
        ax.set_xticks(agg_index[regions[0]:regions[1]])
        ax.set_xticklabels(aggregated["xlabel"][regions[0]:regions[1]])
        for i in range(len(aggregated)):
            if aggregated["color"][aggregated.index[i]] == 'grey':
                ax.bar(agg_index[i],-0.1, yerr = aggregated["stdev"][aggregated.index[i]], capsize=2, width=0.8,edgecolor='black', alpha = 0.5, color=aggregated["color"][aggregated.index[i]], align='center')
            else:
                ax.bar(agg_index[i], aggregated["mean"][aggregated.index[i]], yerr = aggregated["stdev"][aggregated.index[i]], edgecolor='grey',capsize=2, width=0.8, color=aggregated["color"][aggregated.index[i]], align='center')
        ax.axhline(y=0.0, color="silver", linestyle="-")
        legend_elements = [
        Patch(facecolor=ReactivityThreshold.COLOR_HIGH,edgecolor='grey', label="High Reactivity"),
        Patch(facecolor=ReactivityThreshold.COLOR_MEDIUM, edgecolor='grey',label="Medium Reactivity"),
        Patch(facecolor=ReactivityThreshold.COLOR_LOW, edgecolor='grey',label="Low Reactivity"),
        Patch(facecolor=ReactivityThreshold.COLOR_INVALID, alpha = 0.5, edgecolor='black',label="Undetermined"),
        ]
        ax.legend(handles=legend_elements, loc="upper left", ncol=16)
        ax.set_title(f"Nucleotides {aggregated['xlabel'][aggregated.index[regions[0]]].split()[0]} - {aggregated['xlabel'][aggregated.index[regions[1]-1]].split()[0]}",loc='left', fontsize=20)
        ax.set_xlabel("Sequence", fontsize=15)
        ax.set_ylabel("Average reactivity", fontsize=15)
        ax.set_xlim([agg_index[regions[0]] - 1, agg_index[regions[1]-1] + 1])
        ax.set_ylim([min(np.nanmin(aggregated["mean"]-aggregated["stdev"])*1.1,-0.15), np.nanmax(aggregated["mean"]+aggregated["stdev"])*1.2])
        ax2 = ax.twinx()
        ax2.set_ylim([min(np.nanmin(aggregated["mean"]-aggregated["stdev"])*1.1,-0.15), np.nanmax(aggregated["mean"]+aggregated["stdev"])*1.2])
    plt.margins(0)
    plt.suptitle(title, fontsize=30)
    plt.tight_layout()

    
    try:
        plt.tight_layout()
        plt.savefig(output, format=format)
    except ValueError as e:
        logging.error(f"Unable to save plot: {e}")
        open(output, "a").close()

    # mplcursors.cursor()
    # return fig, axs
    # ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right');


def check_files(src, dest):
    if os.path.exists(dest) and os.path.isdir(dest):
        raise fire.core.FireError(
            "Output {0} is a directory ,choose a filename".format(dest)
        )

    for file in src:
        if not os.path.exists(file):
            raise fire.core.FireError('Input file "{0}" does not exists'.format(file))

    return (src, dest)


class ShapeReactivitySeq:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.name = os.path.splitext(os.path.basename(filepath))[0]
        self.df = pd.read_csv(filepath, sep="\t")
        self.df = self.df.set_index("seqNum")
        self.df = self.df.set_index("seqRNA", append=True)
        self.df = self.df.rename_axis(index={"seqRNA": "sequence"})


#def min_enough_values(nvalues: int, min_nsubdata_perc: float = 0.66):
#    if nvalues < 2:
#        return 2
#    return np.ceil(nvalues * min_nsubdata_perc)


#def dispersion_threshold(
#    curmean, max_mean_perc: float = 0.682, min_dispersion: float = 0.5
#):
#    return (
#        curmean * max_mean_perc
#        if curmean * max_mean_perc > min_dispersion
#        else min_dispersion
#    )


def most_uniform_subsample_mean_std(sample):
    smean = np.Infinity
    sstdev = np.Infinity
    ssem = np.Infinity
    smad = np.Infinity
    ssample = None
    for curssample in itertools.combinations(sample, len(sample) - 1):
        curssample = pd.Series(curssample)
        curstdev = curssample.std(ddof=ddof)
        if curstdev < sstdev:
            sstdev = curstdev
            ssem = curssample.sem(ddof=ddof)
            smean = np.mean(curssample)
            smad = (
                (curssample - curssample.mean()).abs().mean()
            )  # Deprecated : curssample.mad()
            ssample = curssample
    return (ssample, smean, sstdev, ssem, smad)


# Alternative method, not used
#def _aggregate_replicates_tscore(row, max_dispersion=0.3):
#    mean = np.NaN
#    stdev = np.NaN
#    desc = "non-consistant"
#    values = row.drop("nrep_with_value").replace(-10, np.NaN)
#    nvalues = values.count()
#    used_values = 0
#
#    # if there is enough values available, we try to compute mean and stdev
#    if values.count() >= row["nrep_with_value"] / 2:
#        curmean = values.mean()
#        curstdev = values.std(ddof=ddof)
#        if curstdev <= max_dispersion:
#            mean = curmean
#            stdev = curstdev
#            used_values = nvalues
#            desc = "accepted"
#        else:
#            subvals = values
#            tscores = subvals.apply(
#                lambda x: abs((x - curmean) / (curstdev / np.sqrt(nvalues)))
#            )
#            while subvals.count() > min_enough_values(row["nrep_with_value"]):
#                tsmaxidx = tscores.idxmax()
#                tscores = tscores.drop(tsmaxidx)
#                subvals = subvals.drop(tsmaxidx)
#
#                curmean = subvals.mean()
#                curstdev = subvals.std()
#                if curstdev <= max_dispersion:
#                    mean = curmean
#                    stdev = curstdev
#                    used_values = subvals.count()
#                    desc = "reduced"
#                    break
#                tscores = subvals.apply(
#                    lambda x: abs((x - curmean) / (curstdev / np.sqrt(subvals.count())))
#                )
#    else:
#        desc = "no-enough-values"
#
#    return pd.Series(
#        {
#            "mean": mean,
#            "stdev": stdev,
#            "used_values": used_values,
#            "desc": desc,
#        },
#        index=["mean", "stdev", "used_values", "desc"],
#    )
#
#
#def aggregate_replicates(
#    row,
#    max_mean_perc: float = 0.682,
#    min_ndata_perc: float = 0.5,
#    min_nsubdata_perc: float = 0.66,
#    min_dispersion: float = 0.05,
#):
#    mean = np.NaN
#    stdev = np.NaN
#    sem = np.NaN
#    mad = np.NaN
#    values = row.drop("nvalid_values").replace(-10, np.NaN).dropna()
#    nvalues = values.count()
#    used_values = 0
#    desc = "non-consistant"
#
#    # Only one value -- no average possible
#    if row["nvalid_values"] == 1 and nvalues == 1:
#        desc = "one-value-available"
#        mean = values[0]
#        stdev = np.NaN
#        sem = np.NaN
#        mad = np.NaN
#        used_values = nvalues
#    else:
#        # if there is enough values available, we try to compute mean and stdev
#        if nvalues >= row["nvalid_values"] * min_ndata_perc:
#            curmean = values.mean()
#            curstdev = values.std(ddof=ddof)
#            cursem = values.sem(ddof=ddof)
#            curmad = (values - values.mean()).abs().mean()  # Deprecated : values.mad()
#
#            if curstdev <= dispersion_threshold(curmean, max_mean_perc, min_dispersion):
#                mean = curmean
#                stdev = curstdev
#                sem = cursem
#                mad = curmad
#                used_values = nvalues
#                desc = "accepted"
#            else:
#                subsample = values
#                while len(subsample) > min_enough_values(
#                    row["nvalid_values"], min_nsubdata_perc
#                ):
#                    (
#                        subsample,
#                        curmean,
#                        curstdev,
#                        cursem,
#                        curmad,
#                    ) = most_uniform_subsample_mean_std(subsample)
#                    if curstdev <= dispersion_threshold(
#                        curmean, max_mean_perc, min_dispersion
#                    ):
#                        mean = curmean
#                        stdev = curstdev
#                        sem = cursem
#                        mad = curmad
#                        used_values = len(subsample)
#                        desc = "reduced"
#                        break
#                if desc == "non-consistant":
#                    mean = -10
#                    stdev = values.std(ddof=ddof)
#                    sem = values.sem(ddof=ddof)
#                    mad = (values - values.mean()).abs().mean()
#        else:
#            desc = "no-enough-values"
#            mean = -10 if any([v == -10 for v in row]) else np.NaN
#    return pd.Series(
#        {
#            "mean": mean,
#            "stdev": stdev,
#            "sem": sem,
#            "mad": mad,
#            "used_values": used_values,
#            "desc": desc,
#        },
#        index=["mean", "stdev", "sem", "mad", "used_values", "desc"],
#    )


def aggregate_replicates(
    row,
    min_std,
    reactivity_medium,
    reactivity_high
):
    mean = np.NaN
    stdev = np.NaN
    sem = np.NaN
    mad = np.NaN
    # values = row.drop("nvalid_values").replace(-10, np.NaN).dropna()
    values = row.drop("nvalid_values").dropna().apply(lambda x: 0 if -1 <= x < 0 else x).apply(lambda x: -10 if -10 < x < -1 else x)
    nvalues = values.replace(-10, np.NaN).dropna().count()
    used_values = 0
    desc = "non-consistant"

    # Only one value -- no average possible
    if row["nvalid_values"] == 1 and nvalues == 1:
        desc = "one-value-available"
        mean = values.iloc[0]
        stdev = np.NaN
        sem = np.NaN
        mad = np.NaN
        used_values = nvalues
    elif list(values).count(-10) < len(list(values))/2:
        nvalid_values = values.replace(-10, np.NaN).dropna()
        mean = nvalid_values.mean()
        stdev = nvalid_values.std(ddof=ddof)
        sem = nvalid_values.sem(ddof=ddof)
        mad = (nvalid_values - nvalid_values.mean()).abs().mean()  # Deprecated : values.mad()

        mean_values = []
        nvalid_values = list(nvalid_values)
        for v in range(len(nvalid_values)-1):
            for u in range(len(nvalid_values)):
                if u > v:
                    mean_uv = (nvalid_values[v] + nvalid_values[u])/2
                    if 0 <= mean < reactivity_medium:
                        if mean_uv >= reactivity_medium:
                            mean_values.append(mean_uv)
                    elif reactivity_medium <= mean < reactivity_high:
                        if reactivity_medium > mean_uv or mean_uv >= reactivity_high:
                            mean_values.append(mean_uv)
                    elif mean >= reactivity_high:
                        if reactivity_high > mean_uv:
                            mean_values.append(mean_uv)
        if mean_values == [] or stdev <= min_std:
            desc = "accepted"
            used_values = nvalues
        elif all(0 <= m < reactivity_medium for m in mean_values) or all(reactivity_medium <= m < reactivity_high for m in mean_values) or all(reactivity_high <= m for m in mean_values):
            desc = 'warning'
            used_values = nvalues
        else:
            mean = -10
            desc = "non-consistant"
    else:
            mean = -10
            desc = "undetermined"

    return pd.Series(
        {
            "mean": mean,
            "stdev": stdev,
            "sem": sem,
            "mad": mad,
            "used_values": used_values,
            "desc": desc,
        },
        index=["mean", "stdev", "sem", "mad", "used_values", "desc"],
    )

def check_duplicated(aggregated: pd.DataFrame):
    dup = pd.DataFrame(aggregated)
    dup = aggregated.reset_index(drop=False)[["seqNum", "sequence"]]

    duplist = dup[dup.duplicated("seqNum", keep=False)]
    dupcount = len(dup[dup.duplicated("seqNum", keep="first")])
    if dupcount > 0:
        raise ValueError(
            f"{dupcount} base/positions are duplicated in aggregated data:"
            " Did you use the same"
            " sequence in QuShape for each replicate ?\n"
            f"bases : {list(duplist)})\n"
        )


def aggregate(
    *files: [str],
    output: str,
    refseq: str = None,
    shape_output=None,
    shape_IP_output=None,
    map_output=None,
    relation_output=None,
    normcol="simple_norm_reactivity",
    min_std,
    reactivity_medium,
    reactivity_high,
#    min_ndata_perc: float = 0.5,
#    min_nsubdata_perc: float = 0.66,
#    max_mean_perc: float = 0.682,
#    min_dispersion: float = 0.05,
    fullplot: str = None,
    plot: str = None,
    plot_title: str = None,
    err_on_dup: bool = True,
):
    """Aggregate reactivity files together

    Aggregate reactivity between all file in output. Calculate mean and stdev
    for coherent value of reactivity

    Parameters
    ----------
    files : [str]
        Files to aggregate together
    output : str
        Output aggregated file
    refseq : str
        Sequence file name
    ipanemap_output :
        Output aggregated file in a compatible format for RNAFold and IPANEMAP
    normcol :
        name of the normalization column in each input file
    min_ndata_perc : float
        (default: 0.5) minimum percentage
        (usable data (not -10)/available data).
        below this value, data row is discarded
    min_nsubdata_perc : float
        (default: 0.66) When try to find consistant mean on a subsample,
        this value represent minimum (subsample size / avalaible data).
        if below this value, data row is discarded
    max_mean_perc : float
        (default: 0.682) in order to be considered as consistant,
        stdev of a row must be below a certain percentage of mean of the row.
        max_mean_perc represent this mean percentage.
    min_dispersion : float
        (default: 0.05) the consistancy ratio cannot be below a certain value

    """

    src = files
    dest = output
    logging.basicConfig(format=f"{dest}: %(levelname)s:%(message)s")

    check_files(src, dest)

    sequence = None
    if refseq is not None:
        import skbio
        fa = skbio.Sequence.read(refseq, format="fasta")
        sequence = str(fa)

    shape_react_seqs = [ShapeReactivitySeq(filepath) for filepath in src]
    shape_dfs = [
        srs.df[[normcol]].rename(columns={normcol: srs.name})
        for srs in shape_react_seqs
    ]

    # print(shape_dfs[2])
    reacts = pd.concat(shape_dfs, axis=1)
    # Remove leading and trailing empty data

    # reacts = reacts.sort_index()
    # first_idx = reacts.first_valid_index()
    # last_idx = reacts.last_valid_index()
    # reacts = reacts.loc[first_idx:last_idx, :]

    reacts["nvalid_values"] = reacts.count(axis=1)

    aggregated = reacts.copy()

    aggregated[
        ["mean", "stdev", "sem", "mad", "used_values", "desc"]
    ] = aggregated.apply(
        lambda row: aggregate_replicates(
            row,
            min_std,
            reactivity_medium,
            reactivity_high
#            max_mean_perc,
#            min_ndata_perc,
#            min_nsubdata_perc,
#            min_dispersion,
        ),
        axis=1,
        result_type="expand",
    )

    aggregated = aggregated[aggregated.columns.drop(["nvalid_values"])]
    aggregated = aggregated.sort_index()
    try:
        check_duplicated(aggregated)
    except ValueError as err:
        logging.error(err.args[0])
        logging.error(f"inputs: {files}")
        if err_on_dup:
            exit(1)

    # Write aggregate.tsv file to results folder
    aggregated.to_csv(dest, sep="\t", float_format="%.4f")

    # ----------------
    # insert/pad missing rows in aggregated
    #
    aggregated.reset_index(inplace=True)
    aggregated.set_index(['seqNum'], inplace=True)

    ## default values for missing rows (others: NaN)
    undef_row = {
        "mean": -10,
        "stdev": 0,
        "used_values": 0,
        "desc": "undetermined"
    }

    ## fill in missing rows
    aggregated_padded = aggregated.reindex(
        list(range(1,aggregated.index[-1]+1))
        ).fillna(undef_row)

    ## set correct sequence
    for i,row in aggregated_padded.iterrows():
        x = row['sequence']
        if type(x)!=str:
            x='N'
            if sequence is not None and 1<=i<=len(sequence):
                x = sequence[i-1]
        aggregated_padded.at[i,'sequence'] = x

    if shape_output is not None:
        aggripan = aggregated_padded[["mean"]]
        aggripan.to_csv(shape_output, sep="\t", float_format="%.4f", header=False)

    if shape_IP_output is not None:
        IP_output = copy.deepcopy(aggregated_padded)
        IP_output.loc[IP_output['desc'] == 'warning', 'mean'] = -10
        aggripan = IP_output[["mean"]]
        aggripan.to_csv(shape_IP_output, sep="\t", float_format="%.4f", header=False)

    if map_output is not None:
        aggripan = aggregated_padded[["mean", "stdev", "sequence"]]
        aggripan.to_csv(map_output, sep="\t", float_format="%.4f", header=False)

    if relation_output is not None:
        rel_output = copy.deepcopy(aggregated)
        rel_output.set_index(['sequence'],inplace=True, append=True)
        rel_output = rel_output.iloc[:,:rel_output.columns.get_loc('mean')]
        correlation_output = rel_output.corr(method='spearman')
        
        sp_sum_corre = {}
        p_sum_corre = {}
        for replicat in correlation_output:
            sp_sum_corre[replicat] = []
            p_sum_corre[replicat] = []
            for cv in correlation_output[replicat]:
                if cv < 0.7:
                    sp_sum_corre[replicat].append(cv)

        rel_output = rel_output.dropna()
        linear_regression = {}
        for i in range(len(rel_output.columns)):
            for j in range(len(rel_output.columns)):
                if j>i:
                    name = rel_output.columns[j] + '='+ 'f('+ rel_output.columns[i] + ')'
                    try:
                        slope, intercept, r_value, p_value, std_err = linregress(rel_output.iloc[:,i],rel_output.iloc[:,j])
                        linear_regression[name] = {'slope':slope,'intercept':intercept,'r_value':r_value, 'r_squared':r_value**2,'p_value':p_value,'std_err':std_err}
                        if r_value < 0.7:
                            p_sum_corre[rel_output.columns[j]].append(r_value)
                            p_sum_corre[rel_output.columns[i]].append(r_value)

                    except ValueError:
                        print(f"Warning: cannot compute linear regression {name}")

        linear_regression = pd.DataFrame(linear_regression).T
                

        with open(relation_output, 'w') as f:
            f.write("Spearman Correlation Matrice\n")
            correlation_output.to_csv(f, sep="\t", float_format="%.4f")
            f.write("\nLinear Regression\n")
            linear_regression.to_csv(f, sep="\t", float_format="%.4f", mode='a')

        sp_warning_replicate = []
        for r in sp_sum_corre:
            if len(sp_sum_corre[r]) > 0.5 * (len(sp_sum_corre)-1):
                sp_warning_replicate.append(r.split('.')[0])

        p_warning_replicate = []
        for r in p_sum_corre:
            if len(p_sum_corre[r]) > 0.5 * (len(sp_sum_corre)-1):
                p_warning_replicate.append(r.split('.')[0])

        if sp_warning_replicate != []:
            print(f"According to the Spearman correlation\n  Please check the input data for replicate(s) : \n{', '.join(sp_warning_replicate)}")
            with open(relation_output, 'a') as f:
                f.write("\n\n\n")
                f.write(f"According to the Spearman correlation\n  Please check the input data for replicate(s) : \n{', '.join(sp_warning_replicate)}\n\n")

        if p_warning_replicate != []:
            print(f"According to the Pearson correlation\n  Please check the input data for replicate(s) : \n{', '.join(p_warning_replicate)}")
            with open(relation_output, 'a') as f:
                f.write("\n\n\n")
                f.write(f"According to the Pearson correlation\n  Please check the input data for replicate(s) : \n{', '.join(p_warning_replicate)}\n\n")

        # correlation_output.to_csv(relation_output, sep="\t", float_format="%.4f")
        # linear_regression.to_csv(relation_output, mode='a', sep="\t", float_format="%.4f")


    if plot is not None or fullplot is not None:
        aggregated.set_index(['sequence'],inplace=True, append=True)
        try:
            title = plot_title if plot_title is not None else output
            plot_aggregate(aggregated, reactivity_medium, reactivity_high, fulloutput=fullplot, output=plot, title=title)
        except ValueError as e:
            logging.error(f"Unable to save plot: {e}")
            open(fullplot, "a").close()
            open(plot, "a").close()


def main():
    return fire.Fire(aggregate)


if __name__ == "__main__":
    main()