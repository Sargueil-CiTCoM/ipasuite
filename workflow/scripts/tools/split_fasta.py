#!/usr/bin/env python3
import fire
import os
from skbio import DNA, RNA


def check_files(src: str, dest: str):
    if os.path.exists(dest) and os.path.isdir(dest):
        raise fire.core.FireError(
            "Output {0} is a directory ,choose a filename".format(dest)
        )

    if not os.path.exists(src):
        raise fire.core.FireError(
            'Input file "{0}" does not exists'.format(src)
        )


def split_fasta(src: str, dest: str, begin: int = 0, end: int = None):
    """split_fasta.

    Parameters
    ----------
    src : str
        Fasta input
    dest : str
        Destination file
    begin : int
        Beginning of fragment
    end : int
        End of fragment
    """

    check_files(src, dest)

    try:
        fa = DNA.read(src, format="fasta")
    except ValueError:
        fa = RNA.read(src, format="fasta")

    sub = fa[begin:end]
    sub.metadata["description"] += f"from {begin} to {end}"

    sub.write(dest, format="fasta")


if __name__ == "__main__":
    fire.Fire(split_fasta)
