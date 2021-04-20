#!/usr/bin/env python3
import sys, os
import pandas as pd
import fire


def get_start_line(filename: str, lookup:str='INDEX') -> int:
    lastlookup = -1
    with open(filename) as file:
        for num, line in enumerate(file, 0):
            if lookup in line: 
                lastlookup = num;
    if lastlookup <= 0:
        raise Exception("invalid file fromat")
    return lastlookup

def to_tsv(src: str, dest: str, lookup:str='INDEX') -> None:
    start = get_start_line(src, lookup)
    df = pd.read_csv(src, skiprows=start, sep="\t")
    df = df.drop(columns=['INDEX', 'CAP'])
    df.to_csv(dest, sep ='\t', index=False)

def all_to_tsv(*files: [str], lookup:str='INDEX', all:bool=False) -> None:
    """Transform seq8000 file into tsv files. last argument is the destination directory"""

    if len(files) < 2:
        raise fire.core.FireError('Needed at least two arguments : [src] and [dest]')
    dest = files[-1]
    src = files[:-1]

    if all:
        if not os.path.isdir(dest):
            raise fire.core.FireError('Output folder {} does not exists', dest)

        for file in src:
            if not os.path.exists(file):
                raise fire.core.FireError('Input file {} does not exists', file)

        for file in src:
            cur_dest = os.path.join(dest, os.path.splitext(os.path.basename(file))[0] + ".tsv");
            to_tsv(file, cur_dest, lookup)
    else:
        if len(files) != 2:
            raise fire.core.FireError('Expected two arguments : [src] and [dest], to treat multiple file, add --all argument')
        to_tsv(src[0], dest, lookup)


if __name__ == "__main__":
    fire.Fire(all_to_tsv)
