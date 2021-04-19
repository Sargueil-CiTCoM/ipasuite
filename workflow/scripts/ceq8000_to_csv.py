#!/usr/bin/env python3
import sys, os
import pandas as pd
import fire


def get_start_line(filename, lookup='INDEX') -> int:
    lastlookup = -1
    with open(filename) as file:
        for num, line in enumerate(file, 0):
            if lookup in line: 
                lastlookup = num;
    if lastlookup <= 0:
        raise Exception("invalid file fromat")
    return lastlookup

def to_csv(src, dest, lookup='INDEX'):
    start = get_start_line(src, lookup)
    df = pd.read_csv(src, skiprows=start, sep="\t")
    df = df.drop(columns=['INDEX', 'CAP'])
    df.to_csv(dest, sep ='\t', index=False)


if __name__ == "__main__":
    fire.Fire(to_csv)
