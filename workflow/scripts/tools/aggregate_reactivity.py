#!/usr/bin/env python3
import sys
import os
import pandas as pd
import fire
import numpy as np

import itertools
from skbio import DNA, RNA
from skbio.alignment import local_pairwise_align_ssw

ddof = 1


def check_files(src, dest):
    if os.path.exists(dest) and os.path.isdir(dest):
        raise fire.core.FireError(
            "Output {0} is a directory ,choose a filename".format(dest)
        )

    for file in src:
        if not os.path.exists(file):
            raise fire.core.FireError(
                'Input file "{0}" does not exists'.format(file)
            )

    return (src, dest)


def get_src_dest(files):
    if len(files) < 2:
        raise fire.core.FireError(
            "Needed at least two arguments : [src] and [dest]"
        )
    dest = files[-1]
    src = files[:-1]

    for file in src:
        if not os.path.exists(file):
            raise fire.core.FireError("Input file {0} does not exists", file)
    if os.path.exists(dest) and os.path.isdir(dest):
        raise fire.core.FireError(
            "Output {0} is a directory ,choose a filename", dest
        )

    return (src, dest)


class ShapeReactivitySeq:
    def __init__(
        self, filepath: str, reference: pd.DataFrame = None, refseq: RNA = None
    ):
        self.filepath = filepath
        self.name = os.path.splitext(os.path.basename(filepath))[0]
        self.df = pd.read_csv(filepath, sep="\t")
        self.sequence = RNA("".join(self.df["seqRNA"]))
        if reference is not None and refseq is not None:
            self.reference = reference
            self.reference_sequence = refseq
            self.init_absolute_position()
        else:
            self.df = self.df.set_index("seqNum")
            self.df = self.df.set_index("seqRNA", append=True)
            self.df = self.df.rename_axis(index={"seqRNA": "sequence"})

    def init_absolute_position(self):
        _, _, startend = local_pairwise_align_ssw(
            self.reference_sequence, self.sequence
        )
        shift = startend[0][0] + 1
        self.df = self.df.set_index(
            [pd.Index(range(shift, len(self.df) + shift), name="abs_position")]
        )


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
    ssample = None
    for curssample in itertools.combinations(sample, len(sample) - 1):
        curssample = pd.Series(curssample)
        curstdev = curssample.std(ddof=ddof)
        if curstdev < sstdev:
            sstdev = curstdev
            ssem = curssample.sem(ddof=ddof)
            smean = np.mean(curssample)
            ssample = curssample
    return (ssample, smean, sstdev, ssem)


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
                    lambda x: abs(
                        (x - curmean) / (curstdev / np.sqrt(subvals.count()))
                    )
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
        used_values = nvalues
    else:
        # if there is enough values available, we try to compute mean and stdev
        if nvalues >= row["nvalid_values"] * min_ndata_perc:
            curmean = values.mean()
            curstdev = values.std(ddof=ddof)
            cursem = values.sem(ddof=ddof)

            if curstdev <= dispersion_threshold(
                curmean, max_mean_perc, min_dispersion
            ):
                mean = curmean
                stdev = curstdev
                sem = cursem
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
                    ) = most_uniform_subsample_mean_std(subsample)
                    if curstdev <= dispersion_threshold(
                        curmean, max_mean_perc, min_dispersion
                    ):
                        mean = curmean
                        stdev = curstdev
                        sem = cursem
                        used_values = len(subsample)
                        desc = "reduced"
                        break
                if desc == "non-consistant":
                    mean = -10
                    stdev = values.std(ddof=ddof)
                    sem = values.sem(ddof=ddof)
        else:
            desc = "no-enough-values"
            mean = -10 if any([v == -10 for v in row]) else np.NaN
    return pd.Series(
        {
            "mean": mean,
            "stdev": stdev,
            "sem": sem,
            "used_values": used_values,
            "desc": desc,
        },
        index=["mean", "stdev", "sem", "used_values", "desc"],
    )


def aggregate(
    *files: [str],
    output: str,
    refseq=None,
    ipanemap_output=None,
    ref_is_dna=False,
    normcol="simple_norm_reactivity",
    min_ndata_perc: float = 0.5,
    min_nsubdata_perc: float = 0.66,
    max_mean_perc: float = 0.682,
    min_dispersion: float = 0.05
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
    refseq :
        file path to the reference sequence, None by default
    ipanemap_output :
        Output aggregated file in a compatible format for RNAFold and IPANEMAP
    ref_is_dna :
        if reference file is DNA, put --ref_is_dna=True, rna by default
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

    check_files(src, dest)

    refseqdf = None
    if refseq:
        if ref_is_dna:
            refseq = DNA.read(refseq, format="fasta").transcribe()
        else:
            try:
                refseq = RNA.read(refseq, format="fasta")
            except ValueError:
                sys.stderr.write(
                    "WARNING: File is not RNA, using DNA + transcription\n"
                )
                refseq = DNA.read(refseq, format="fasta").transcribe()

        refseqdf = pd.DataFrame(
            [v.decode("utf-8") for v in refseq.values], columns=["sequence"]
        )
        refseqdf.index += 1
        refseqdf.index.name = "abs_position"

    shape_react_seqs = [
        ShapeReactivitySeq(filepath, refseqdf, refseq) for filepath in src
    ]
    shape_dfs = []
    if refseqdf is not None:
        shape_dfs = [refseqdf]
        # print(refseqdf)

    shape_dfs.extend(
        [
            srs.df[[normcol]].rename(columns={normcol: srs.name})
            for srs in shape_react_seqs
        ]
    )

    # print(shape_dfs[2])
    reacts = pd.concat(shape_dfs, axis=1)
    if refseqdf is not None:
        reacts = reacts.set_index(["sequence"], append=True)
    # Remove leading and trailing empty data

    # reacts = reacts.sort_index()
    # first_idx = reacts.first_valid_index()
    # last_idx = reacts.last_valid_index()
    # reacts = reacts.loc[first_idx:last_idx, :]

    reacts["nvalid_values"] = reacts.count(axis=1)

    aggregated = reacts.copy()

    aggregated[
        ["mean", "stdev", "sem", "used_values", "desc"]
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

    # aggregated = aggregated.reset_index(level="seqRNA")

    # print(aggregated)
    aggregated.to_csv(dest, sep="\t", float_format="%.4f")
    if ipanemap_output is not None:
        aggregated.reset_index(level="sequence")["mean"].to_csv(
            ipanemap_output, sep="\t", float_format="%.4f", header=False
        )


if __name__ == "__main__":
    fire.Fire(aggregate)
