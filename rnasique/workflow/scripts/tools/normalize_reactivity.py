#!/usr/bin/env python3
import sys
import os
import pandas as pd
import fire
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_style("dark")

allowed_methods = ["simple", "interquartile"]


def plot_norm_reactivity(
    df: pd.DataFrame,
    title="Normalized reactivity",
    output="fig.svg",
    format="svg",
):
    df = df.sort_values(by=["seqNum"], ascending=True)
    df["xlabel"] = df["seqRNA"].astype(str) + "\n" + df.index.astype(str)
    df.loc[df["simple_norm_reactivity"] == -10, "simple_norm_reactivity"] = 0
    df.loc[
        df["interquartile_norm_reactivity"] == -10,
        "interquartile_norm_reactivity",
    ] = 0

    ax = df.plot(
        x="xlabel",
        y=["simple_norm_reactivity", "interquartile_norm_reactivity"],
        # width=0.7,
        rot=0,
        # kind="line",
        drawstyle="steps",
        figsize=(len(df) / 3.5, 4),
        x_compat=True,
        xticks=np.arange(0, len(df) + 1, 1)
    )
    ax.set_xlabel("Sequence")
    ax.set_ylabel("Normalized reactivity")
    plt.margins(0)
    plt.title(title, loc="left")
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(output, format=format)


def check_files(src, dest):
    if len(src) == 0:
        raise fire.core.FireError("missing input")

    if len(src) > 1:
        if dest is None:
            raise fire.core.FireError(
                "--output is mandatory when more than on file to normalize"
            )
        elif not os.path.isdir(dest):
            raise fire.core.FireError(
                'Output folder "{0}" does not exists '
                "or is not a directory".format(dest)
            )
    for file in src:
        if not os.path.exists(file):
            raise fire.core.FireError(
                'Input file "{0}" does not exists'.format(file)
            )

    return (src, dest)


def check_inputs_columns(df: pd.DataFrame) -> None:
    try:
        assert df.dtypes["areaBG"] == np.float64
        assert df.dtypes["areaRX"] == np.float64
        assert df.dtypes["seqNum"] == np.int64
        assert df.dtypes["seqRNA"] == object
    except:
        raise Exception(
            "Invalid input columns:\n"
            "needed at least: {{'areaBG': float64,"
            "'areaRX': float64, 'seqNum': np.int64, 'seqRNA: np.string'}}\n"
            "got: {}".format(df.dtypes)
        )


def check_norm_methods(methods: [str]) -> None:
    if not all((method in allowed_methods) for method in methods):
        raise Exception(
            "invalid normalization method : allowed methods are {0}",
            allowed_methods,
        )


def prepare_normalization_df(
    df: pd.DataFrame, reactive_nucleotides: [str] = ["A", "C", "G", "U"]
) -> pd.DataFrame:

    # Keep only reactive nucleotides for
    df = df.loc[
        (
            (df["seqRNA"].isin(reactive_nucleotides))
            & (df["corr_areaRX"].notna())
        ),
        :,
    ].copy(deep=True)
    return df


def compute_simple_norm_term(
    df: pd.DataFrame,
    norm_column: str = "corr_areaRX",
    stop_percentile: float = 90,
    outlier_percentile: float = 98,
    norm_term_avg_percentile: float = 90,
) -> float:
    """\"Simple\" normalization

    average of the 10% top reactive nucleotides, minus the top 2% (outliers)

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe used to calculate normalization term
    norm_column : str
       df column used to calculate normalization term
    stop_percentile : float
         The threshold above which background is estimated to be too high -
         data above this threshold are excluded from nomalization
    outlier_percentile : float
         threshold (in percent) above which reactivity
         is considered as too high
    norm_term_avg_percentile : float
        threshold (in percent) above which reactivities are
        used as to calculate normalization term

    Returns
    -------
    float

    """

    # Get value between perc(norm_term_avg:90) perc(outlier:98)
    df = df.loc[
        (
            df["corr_areaRX"]
            > np.percentile(df["corr_areaRX"], norm_term_avg_percentile)
        )
        & (
            df["corr_areaRX"]
            < np.percentile(df["corr_areaRX"], outlier_percentile)
        )
    ]

    # Average this value
    norm_term = np.average(df["corr_areaRX"])
    return norm_term


def compute_interquart_norm_term(
    df: pd.DataFrame, norm_column: str = "corr_areaRX"
) -> float:
    """\"BoxPlot\"  normalization

        define outline threshold as 1.5 * interquartile range
        remove value above 3rd quartile that are above threshold,
        average arround 10% most reactive position remaining

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe used to calculate normalization term
    norm_column : str
        column name used for normalization

    Returns
    -------
    float
        normalization term for df

    """

    # Get threshold outliers as 1.5*(val(75%) - val(25%))
    interquartile_range = np.percentile(df[norm_column], 75) - np.percentile(
        df[norm_column], 25
    )

    # Normprobing filtering

    interq_thres = 1.5 * interquartile_range
    # Weeks 2010 filtering
    # Removing value above  3rd quartile
    # + interquartile_range * 1.5 -- Q3 + 1.5 (Q3 - Q1)
    df_no_outlier = df.loc[
        (df[norm_column] < np.percentile(df[norm_column], 75) + interq_thres)
    ]

    # Norm probing filtering
    # Removing value above threshold and value above 3rd quartile
    # equiv to : keep value below threshold AND below 3rd quartile

    # df_no_outlier = df.loc[(df[norm_column] < interq_thres) &
    # (df[norm_column] < np.percentile(df[norm_column], 75))]

    most_reactive_df = df_no_outlier.loc[
        df_no_outlier[norm_column]
        > np.percentile(df_no_outlier[norm_column], 90)
    ]
    norm_term = np.average(most_reactive_df[norm_column])
    return norm_term


def normalize_one_path(
    inputpath: str,
    outputpath: str = None,
    reactive_nucleotides: [str] = ["G", "C", "A", "U"],
    stop_percentile: float = 90.0,
    simple_outlier_percentile: float = 98.0,
    simple_norm_term_avg_percentile: float = 90.0,
    low_norm_reactivity_threshold: float = -0.3,
    norm_methods: [str] = ["simple", "interquartile"],
    normcol: str = "simple_norm_reactivity",
    plot: str = None,
    plot_title: str = None,
    shape_output: bool = None,
    map_output: bool = None,

) -> int:

    intensity_area_df = pd.read_csv(inputpath, sep="\t")
    check_inputs_columns(intensity_area_df)

    # Adding seqNum as index
    intensity_area_df = (
        intensity_area_df[["seqNum", "seqRNA", "areaRX", "areaBG"]]
        .set_index("seqNum", drop=True)
        .sort_index()
    )

    # To 10% of areaBG are considered as bad quality data,
    # and will not be used.

    # background correction
    intensity_area_df.loc[:, "corr_areaRX"] = (
        intensity_area_df["areaRX"] - intensity_area_df["areaBG"]
    )

    intensity_area_df.loc[
        intensity_area_df["areaBG"]
        > np.percentile(intensity_area_df["areaBG"], stop_percentile),
        "corr_areaRX",
    ] = np.NaN

    norm_df = intensity_area_df.copy(deep=True)

    if "simple" in norm_methods:
        df = prepare_normalization_df(norm_df)
        simple_norm_term = compute_simple_norm_term(
            df,
            outlier_percentile=simple_outlier_percentile,
            norm_term_avg_percentile=simple_norm_term_avg_percentile,
        )
        df["simple_norm_reactivity"] = df["corr_areaRX"] / simple_norm_term
        df.loc[
            (df["simple_norm_reactivity"] < low_norm_reactivity_threshold),
            "simple_norm_reactivity",
        ] = 0

        norm_df["simple_norm_reactivity"] = np.NaN
        # norm_df.loc[norm_df['seqRNA'].isin(reactive_nucleotides),
        # "simple_norm_reactivity"] = 0
        norm_df.loc[df.index, "simple_norm_reactivity"] = df[
            "simple_norm_reactivity"
        ]

    if "interquartile" in norm_methods:
        df = prepare_normalization_df(norm_df)
        inter_norm_term = compute_interquart_norm_term(df)
        df["interquartile_norm_reactivity"] = (
            df["corr_areaRX"] / inter_norm_term
        )
        df.loc[
            (
                df["interquartile_norm_reactivity"]
                < low_norm_reactivity_threshold
            ),
            "interquartile_norm_reactivity",
        ] = 0

        norm_df["interquartile_norm_reactivity"] = np.NaN
        # norm_df.loc[norm_df['seqRNA'].isin(reactive_nucleotides),
        # "interquartile_norm_reactivity"] = 0
        norm_df.loc[df.index, "interquartile_norm_reactivity"] = df[
            "interquartile_norm_reactivity"
        ]

    norm_df.fillna(-10, inplace=True)
    if outputpath is None:
        outputpath = sys.stdout
    norm_df.to_csv(outputpath, sep="\t", float_format="%.4f")

    if shape_output is not None:
        normipan = norm_df[[normcol]]
        idxmin = normipan.index.min()
        firstrows = pd.DataFrame(
            {normcol: np.full(idxmin - 1, -10)}, index=range(1, idxmin)
        )
        firstrows.index.names = ["seqNum"]
        normipan = pd.concat([firstrows, normipan])
        normipan.to_csv(shape_output, sep="\t", float_format="%.4f", header=False)

    if map_output is not None:
        norm_df["stdev"] = 0
        normipan = norm_df[[normcol, "stdev", "seqRNA"]]
        idxmin = normipan.index.min()
        firstrows = pd.DataFrame(
            {
                normcol: np.full(idxmin - 1, -10),
                "stdev": np.zeros(idxmin - 1),
                "seqRNA": np.full(idxmin - 1, "N"),
            },
            index=range(1, idxmin),
        )
        firstrows.index.names = ["seqNum"]
        normipan = pd.concat([firstrows, normipan])
        normipan["stdev"] = 0
        normipan.to_csv(map_output, sep="\t", float_format="%.4f", header=False)

    if plot is not None:
        title = plot_title if plot_title is not None else outputpath
        plot_norm_reactivity(norm_df, title, plot)

    return 0


def normalize_all(
    *inputs: [str],
    output: str = None,
    reactive_nucleotides: [str] = ["G", "C", "A", "U"],
    stop_percentile: float = 90.0,
    simple_outlier_percentile: float = 98.0,
    simple_norm_term_avg_percentile: float = 90.0,
    low_norm_reactivity_threshold: float = -0.3,
    norm_methods: [str] = ["simple", "interquartile"],
    normcol: str = "simple_norm_reactivity",
    plot: str = None,
    plot_title: str = None,
    shape_output: bool = None,
    map_output: bool = None,
) -> int:
    """Normalized reactivity for each input files

       Output tsv file with normalized reactivities

    Parameters
    ----------
    inputs : [str]
        List of file to normalize in tsv
    output : str
        Path directory where to output .tsv files
        containing normalized reactivity
    reactive_nucleotides : [str]
        (default: A,C,G,U) comma separeted list of reactive nucleotides
        (A,C,G,U) are accepted
    stop_percentile : float
        (default: 90. )The threshold above which background is estimated
        to be too high - data above this threshold will be discarded
    simple_outlier_percentile : float
        (default)simple method only - threshold (in percent) above which
        reactivity is considered as too high
    simple_norm_term_avg_percentile : float
        simple_method_only - threshold (in percent) above which reactivities
        are used as to calculate normalization term
    low_norm_reactivity_threshold : float
        normalized reactivity threshold above which reactivity
        is not considered as significant, and removed
    norm_methods : [str]
        comma-separated list of normalization methods:
        `simple` and `interquartile` are allowed

    Returns
    -------
    int
        output tsv files

    """
    res = 0
    check_norm_methods(norm_methods)
    check_files(inputs, output)

    kwargs = {
        "reactive_nucleotides": reactive_nucleotides,
        "stop_percentile": stop_percentile,
        "simple_outlier_percentile": simple_outlier_percentile,
        "simple_norm_term_avg_percentile": simple_norm_term_avg_percentile,
        "low_norm_reactivity_threshold": low_norm_reactivity_threshold,
        "norm_methods": norm_methods,
        "normcol": normcol,
        "plot": plot,
        "plot_title": plot_title,
        "shape_output": shape_output,
        "map_output": map_output
    }

    if output is None or not os.path.isdir(output):
        res += normalize_one_path(inputs[0], output, **kwargs)
    else:
        for file in inputs:
            cur_dest = os.path.join(
                output,
                os.path.splitext(os.path.basename(file))[0] + ".tsv.txt",
            )
            res += normalize_one_path(output, cur_dest, **kwargs)

    # return res


def main():
    return fire.Fire(normalize_all)


if __name__ == "__main__":
    main()
