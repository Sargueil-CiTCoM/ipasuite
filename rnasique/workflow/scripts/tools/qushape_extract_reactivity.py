#!/usr/bin/env python3

from .utils import qushapeFuncReport as qsfr
from .utils import fasta
import os
import fire
import bsddb3 as bsddb
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Hack to access data from ill formed .qushape file
sys.path.append(os.path.join(os.path.dirname(__file__), "qushape_helper"))

plt.style.library["seaborn-dark"]


def plot_reactivity(
    df: pd.DataFrame, title="Reactivity", output="fig.svg", format="svg"
):
    df = df.sort_values(by=["seqNum"], ascending=True)
    df["xlabel"] = df["seqRNA"].astype(str) + "\n" + df["seqNum"].astype(str)

    ax = df.plot(
        x="xlabel",
        y=["areaRX", "areaBG", "areaDiff"],
        width=0.7,
        rot=0,
        kind="bar",
        figsize=(len(df), 4),
    )
    ax.set_xlabel("Sequence")
    ax.set_ylabel("Raw reactivity")
    plt.title(title, loc="left")
    plt.legend(loc="upper left")

    plt.tight_layout()
    plt.savefig(output, format=format)


def getProjData(filepath):
    proj = None
    if os.path.exists(filepath):

        db = bsddb.hashopen(filepath)
        try:
            proj = pickle.loads(db[b"dProject"], encoding="latin1")
        except UnicodeDecodeError:
            proj = pickle.loads(db[b"dProject"])

        db.close()
    else:
        raise ValueError("Project file '{0}' not found".format(filepath))
    return proj


def check_qushape_using_correct_rna(project, rna_file):

    seqname, seq = fasta.get_first_fasta_seq(rna_file)
    seq = seq.replace("T", "U")
    if not seq.startswith(project["RNA"]):
        print(
            "Qushape project RNA is different from provided rna_file :\n"
            f"qushape seq ({project['fNameSeq']}) : {project['RNA']}\n"
            f"rna_file seq {seqname} - {rna_file}: {seq}"
        )
        raise SystemExit(-1)


def extract_reactivity(
    qushape_project: str,
    output: str = "reactivity.tsv.txt",
    plot: str = None,
    plot_title: str = None,
    rna_file: str = None,
    dry_run=False,
):

    proj = getProjData(qushape_project)

    if rna_file is not None:
        check_qushape_using_correct_rna(proj, rna_file)

    #    if not "Reactivity" in proj["scriptList"][-1]:
    if "area" not in proj["dPeakRX"] or len(proj["seqNum"]) == 0:
        print("You must finish QuShape treatment to extract reactivity")
        raise SystemExit(-1)

    report = qsfr.createDReport(proj)
    if not dry_run:
        qsfr.writeReportFile(report, output)

    if plot is not None:
        df = pd.read_csv(output, sep="\t")
        title = plot_title if plot_title is not None else output
        plot_reactivity(df, title, plot)


def main():
    return fire.Fire(extract_reactivity)


if __name__ == "__main__":
    main()
