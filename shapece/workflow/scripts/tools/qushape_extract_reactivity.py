#!/usr/bin/env python3
import utils.qushapeFuncReport as qsfr
import os
import fire
import bsddb3 as bsddb
import pickle
import pandas as pd
import matplotlib.pyplot as plt


def plot_reactivity(
    df: pd.DataFrame, title="Reactivity", output="fig.svg", format="svg"
):
    df = df.sort_values(by=["seqNum"], ascending=True)
    df["xlabel"] = df["seqRNA"].astype(str) + "\n" + df["seqNum"].astype(str)

    df.plot(
        x="xlabel",
        y=["areaRX", "areaBG", "areaDiff"],
        width=0.7,
        rot=0,
        kind="bar",
        figsize=(len(df), 4),
    )
    plt.title(title, loc='left')
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig(output, format=format)


def getProjData(filepath):
    proj = None
    if os.path.exists(filepath):

        db = bsddb.hashopen(filepath)
        try:
            proj = pickle.loads(db[b"dProject"], encoding="latin1")
        except:
            proj = pickle.loads(db[b"dProject"])

        db.close()
    else:
        raise ValueError("Project file '{0}' not found".format(filepath))
    return proj


def extract_reactivity(
    qushape_project: str,
    output: str = "reactivity.tsv.txt",
    plot: str = None,
    dry_run=False,
):

    proj = getProjData(qushape_project)

    #    if not "Reactivity" in proj["scriptList"][-1]:
    if "area" not in proj["dPeakRX"] or len(proj["seqNum"]) == 0:
        print("You must finish QuShape treatment to extract reactivity")
        raise SystemExit(-1)

    report = qsfr.createDReport(proj)
    if not dry_run:
        qsfr.writeReportFile(report, output)

    if plot is not None:
        df = pd.read_csv(output, sep='\t')
        plot_reactivity(df, output, plot)


if __name__ == "__main__":
    fire.Fire(extract_reactivity)
