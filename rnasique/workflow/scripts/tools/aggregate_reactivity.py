#!/usr/bin/env python3
import os
import pandas as pd
import fire
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import cm
import itertools
import logging

plt.style.library["seaborn-dark"]

ddof = 1


class ReactivityThreshold:
    INVALID = -0.3

    LOW = 0.5
    MEDIUM = 0.8
    HIGH = 1.0

    COLOR_INVALID = "grey"
    COLOR_NONE = "white"
    COLOR_LOW = "yellow"
    COLOR_MEDIUM = "orange"
    COLOR_HIGH = "red"


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
        if row["desc"] == "reduced":
            ax.axvspan(
                (idx - first_idx) - 0.5,
                (idx - first_idx + 0.4),
                alpha=0.3,
                color="orange",
                label="_" * fr + "Reduced value",
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
        if row["desc"] == "no-enough-values":
            ax.axvspan(
                (idx - first_idx) - 0.5,
                (idx - first_idx + 0.4),
                alpha=0.3,
                color="blue",
                label="_" * fnev + "Not enough value",
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
        colormap=cm.cubehelix,
        linewidth=0.5,
    )

    plot_aggregation_infos(aggregated, ax)

    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.margins(0)
    plt.title(title, loc="left")
    plt.legend(loc="upper left")
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
    except Exception as e:
        logging.error("fail to generate error bar")
        # logging.error(e)

    ax.set_xlabel("Sequence")
    ax.set_ylabel("Aggregated reactivity")
    try:
        plt.tight_layout()
        plt.savefig(fulloutput, format=format)
    except ValueError as e:
        logging.error(f"Unable to save fullplot: {e}")
        open(fulloutput, "a").close()

    fig, ax = plt.subplots(
        nrows=1, ncols=1, sharex=True, figsize=(len(aggregated) / 4, 4)
    )

    aggregated["color"] = ReactivityThreshold.COLOR_NONE
    aggregated.loc[
        (aggregated["mean"] > ReactivityThreshold.HIGH), "color"
    ] = ReactivityThreshold.COLOR_HIGH
    aggregated.loc[
        (
            (aggregated["mean"] <= ReactivityThreshold.HIGH)
            & (aggregated["mean"] > ReactivityThreshold.MEDIUM)
        ),
        "color",
    ] = ReactivityThreshold.COLOR_MEDIUM
    aggregated.loc[
        (
            (aggregated["mean"] <= ReactivityThreshold.MEDIUM)
            & (aggregated["mean"] > ReactivityThreshold.LOW)
        ),
        "color",
    ] = ReactivityThreshold.COLOR_LOW
    aggregated.loc[
        (aggregated["mean"] < ReactivityThreshold.INVALID), "color"
    ] = ReactivityThreshold.COLOR_INVALID

    aggregated.loc[(aggregated["mean"] == -10), "stdev"] = np.NaN
    aggregated.loc[(aggregated["mean"] == -10), "mean"] = np.NaN
    aggregated["xlabel_rot"] = (
        aggregated.index.get_level_values("seqNum").astype(str)
        + " - "
        + aggregated.index.get_level_values("sequence").astype(str)
    )

    aggregated.plot(
        ax=ax,
        x="xlabel_rot",
        rot=70,
        y="mean",
        kind="bar",
        width=1,
        color=aggregated["color"],
        yerr=None if aggregated["stdev"].isnull().all() else "stdev",
        capsize=2,
        stacked=False,
    )
    ax = meanstdev[["mean", "xlabel"]].plot(
        x="xlabel",
        y="mean",
        drawstyle="steps-mid",
        ax=ax,
        colormap=cm.cubehelix,
        linewidth=0.5,
    )
    ax.set_xlabel("Sequence")
    ax.set_ylabel("Aggregated reactivity")
    plt.margins(0)
    plt.title(title, loc="left")

    legend_elements = [
        Patch(facecolor=ReactivityThreshold.COLOR_HIGH, label="High Reactivity"),
        Patch(facecolor=ReactivityThreshold.COLOR_MEDIUM, label="Medium Reactivity"),
        Patch(facecolor=ReactivityThreshold.COLOR_LOW, label="Low Reactivity"),
        Patch(facecolor=ReactivityThreshold.COLOR_INVALID, label="Invalid Reactivity"),
    ]
    plt.legend(handles=legend_elements, loc="upper left")
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


def min_enough_values(nvalues: int, min_nsubdata_perc: float = 0.66):
    if nvalues < 2:
        return 2
    return np.ceil(nvalues * min_nsubdata_perc)


def dispersion_threshold(
    curmean, max_mean_perc: float = 0.682, min_dispersion: float = 0.5
):
    return (
        curmean * max_mean_perc
        if curmean * max_mean_perc > min_dispersion
        else min_dispersion
    )


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
def _aggregate_replicates_tscore(row, max_dispersion=0.3):
    mean = np.NaN
    stdev = np.NaN
    desc = "non-consistant"
    values = row.drop("nrep_with_value").replace(-10, np.NaN)
    nvalues = values.count()
    used_values = 0

    # if there is enough values available, we try to compute mean and stdev
    if values.count() >= row["nrep_with_value"] / 2:
        curmean = values.mean()
        curstdev = values.std(ddof=ddof)
        if curstdev <= max_dispersion:
            mean = curmean
            stdev = curstdev
            used_values = nvalues
            desc = "accepted"
        else:
            subvals = values
            tscores = subvals.apply(
                lambda x: abs((x - curmean) / (curstdev / np.sqrt(nvalues)))
            )
            while subvals.count() > min_enough_values(row["nrep_with_value"]):
                tsmaxidx = tscores.idxmax()
                tscores = tscores.drop(tsmaxidx)
                subvals = subvals.drop(tsmaxidx)

                curmean = subvals.mean()
                curstdev = subvals.std()
                if curstdev <= max_dispersion:
                    mean = curmean
                    stdev = curstdev
                    used_values = subvals.count()
                    desc = "reduced"
                    break
                tscores = subvals.apply(
                    lambda x: abs((x - curmean) / (curstdev / np.sqrt(subvals.count())))
                )
    else:
        desc = "no-enough-values"

    return pd.Series(
        {
            "mean": mean,
            "stdev": stdev,
            "used_values": used_values,
            "desc": desc,
        },
        index=["mean", "stdev", "used_values", "desc"],
    )


def aggregate_replicates(
    row,
    max_mean_perc: float = 0.682,
    min_ndata_perc: float = 0.5,
    min_nsubdata_perc: float = 0.66,
    min_dispersion: float = 0.05,
):
    mean = np.NaN
    stdev = np.NaN
    sem = np.NaN
    mad = np.NaN
    values = row.drop("nvalid_values").replace(-10, np.NaN).dropna()
    nvalues = values.count()
    used_values = 0
    desc = "non-consistant"

    # Only one value -- no average possible
    if row["nvalid_values"] == 1 and nvalues == 1:
        desc = "one-value-available"
        mean = values[0]
        stdev = np.NaN
        sem = np.NaN
        mad = np.NaN
        used_values = nvalues
    else:
        # if there is enough values available, we try to compute mean and stdev
        if nvalues >= row["nvalid_values"] * min_ndata_perc:
            curmean = values.mean()
            curstdev = values.std(ddof=ddof)
            cursem = values.sem(ddof=ddof)
            curmad = (values - values.mean()).abs().mean()  # Deprecated : values.mad()

            if curstdev <= dispersion_threshold(curmean, max_mean_perc, min_dispersion):
                mean = curmean
                stdev = curstdev
                sem = cursem
                mad = curmad
                used_values = nvalues
                desc = "accepted"
            else:
                subsample = values
                while len(subsample) > min_enough_values(
                    row["nvalid_values"], min_nsubdata_perc
                ):
                    (
                        subsample,
                        curmean,
                        curstdev,
                        cursem,
                        curmad,
                    ) = most_uniform_subsample_mean_std(subsample)
                    if curstdev <= dispersion_threshold(
                        curmean, max_mean_perc, min_dispersion
                    ):
                        mean = curmean
                        stdev = curstdev
                        sem = cursem
                        mad = curmad
                        used_values = len(subsample)
                        desc = "reduced"
                        break
                if desc == "non-consistant":
                    mean = -10
                    stdev = values.std(ddof=ddof)
                    sem = values.sem(ddof=ddof)
                    mad = (values - values.mean()).abs().mean()
        else:
            desc = "no-enough-values"
            mean = -10 if any([v == -10 for v in row]) else np.NaN
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
    shape_output=None,
    map_output=None,
    normcol="simple_norm_reactivity",
    min_ndata_perc: float = 0.5,
    min_nsubdata_perc: float = 0.66,
    max_mean_perc: float = 0.682,
    min_dispersion: float = 0.05,
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

    shape_react_seqs = [ShapeReactivitySeq(filepath) for filepath in src]
    shape_dfs = []

    shape_dfs.extend(
        [
            srs.df[[normcol]].rename(columns={normcol: srs.name})
            for srs in shape_react_seqs
        ]
    )

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
            max_mean_perc,
            min_ndata_perc,
            min_nsubdata_perc,
            min_dispersion,
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
    # aggregated = aggregated.reset_index(level="seqRNA")

    # print(aggregated)
    aggregated.to_csv(dest, sep="\t", float_format="%.4f")
    if shape_output is not None:
        aggripan = aggregated.reset_index(level="sequence")[["mean"]]
        idxmin = aggripan.index.min()
        firstrows = pd.DataFrame(
            {"mean": np.full(idxmin - 1, -10)}, index=range(1, idxmin)
        )
        firstrows.index.names = ["seqNum"]
        aggripan = pd.concat([firstrows, aggripan])
        aggripan.to_csv(shape_output, sep="\t", float_format="%.4f", header=False)

    if map_output is not None:
        aggripan = aggregated.reset_index(level="sequence")[
            ["mean", "stdev", "sequence"]
        ]
        idxmin = aggripan.index.min()
        firstrows = pd.DataFrame(
            {
                "mean": np.full(idxmin - 1, -10),
                "stdev": np.zeros(idxmin - 1),
                "sequence": np.full(idxmin - 1, "N"),
            },
            index=range(1, idxmin),
        )
        firstrows.index.names = ["seqNum"]
        aggripan = pd.concat([firstrows, aggripan])
        aggripan.to_csv(map_output, sep="\t", float_format="%.4f", header=False)
    if plot is not None or fullplot is not None:
        try:
            title = plot_title if plot_title is not None else output
            plot_aggregate(aggregated, fulloutput=fullplot, output=plot, title=title)
        except ValueError as e:
            logging.error(f"Unable to save plot: {e}")
            open(fullplot, "a").close()
            open(plot, "a").close()


def main():
    return fire.Fire(aggregate)


if __name__ == "__main__":
    main()
