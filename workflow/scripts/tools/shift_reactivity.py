#!/usr/bin/env python3
import fire
import os
from skbio import DNA, RNA
import pandas as pd


def check_files(src: str, dest: str):
    if os.path.exists(dest) and os.path.isdir(dest):
        raise fire.core.FireError(
            "Output {0} is a directory ,choose a filename".format(dest)
        )

    if not os.path.exists(src):
        raise fire.core.FireError(
            'Input file "{0}" does not exists'.format(src)
        )


def shift(src: str, reference: str, dest: str, begin: int = 0):
    check_files(src, dest)
    try:
        refseq = DNA.read(reference, format="fasta")
    except ValueError:
        refseq = RNA.read(reference, format="fasta")

    refdf = pd.DataFrame(
        [v.decode("utf-8") for v in refseq.values], columns=["seqRNA"]
    )
    refdf.index += 1
    refdf.index.name = "seqNum"

    srcdf = pd.read_csv(src, sep="\t")
    srcdf = srcdf.rename(columns={"seqNum": "rel_seqNum", "seqRNA": "rel_seqRNA"})
    srcdf["seqNum"] = srcdf["rel_seqNum"] + begin - 1
    srcdf = srcdf.set_index("seqNum")
    print(srcdf)
    print(refdf)
    outdf = pd.concat([refdf, srcdf], axis=1)

    print(outdf)

    outdf.to_csv(dest, sep="\t", float_format="%.4f")


if __name__ == "__main__":
    fire.Fire(shift)
