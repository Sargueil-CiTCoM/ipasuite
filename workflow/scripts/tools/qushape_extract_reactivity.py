#!/usr/bin/env python3
import utils.qushapeFuncReport as qsfr
import os
import fire
import bsddb3 as bsddb
import pickle


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
    qushape_project: str, output: str = "reactivity.tsv.txt", dry_run=False
):

    proj = getProjData(qushape_project)

    #    if not "Reactivity" in proj["scriptList"][-1]:
    if "area" not in proj["dPeakRX"] or len(proj["seqNum"]) == 0:
        print("You must finish QuShape treatment to extract reactivity")
        raise SystemExit(-1)

    report = qsfr.createDReport(proj)
    if not dry_run:
        qsfr.writeReportFile(report, output)


if __name__ == "__main__":
    fire.Fire(extract_reactivity)
