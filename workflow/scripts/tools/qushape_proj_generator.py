#!/usr/bin/python3

import utils.qushapeFuncFile as qsff
import utils.qushapeFuncGeneral as qsfg
import shelve
import bsddb3 as bsddb
import fire
import os
import pickle
from copy import deepcopy


def getProjRefData(filepath):
    refproj = None
    if os.path.exists(filepath):

        db = bsddb.hashopen(filepath)
        try:
            refproj = pickle.loads(db[b"dProject"], encoding="latin1")
        except:
            refproj = pickle.loads(db[b"dProject"])

        db.close()
    else:
        raise ValueError("Ref project file '{0}' not found".format(filepath))
    return refproj


def createDProject(
    rxfile: str,
    bgfile: str,
    output: str = "output.qushape",
    refseq: str = None,
    refproj: str = None,
    channels=None,
    ddNTP="ddC",
    path: str = None,
):

    if path is None:
        path = os.getcwd()
    if channels is None:
        channels = {"RX": 0, "RXS1": 2, "BG": 0, "BGS1": 2}

    if not os.path.exists(rxfile):
        raise ValueError("RX file '{0}' not found".format(rxfile))
    if not os.path.exists(bgfile):
        raise ValueError("BX file '{0}' not found".format(bgfile))
    if refseq is not None and not os.path.exists(refseq):
        raise ValueError("SEQ file '{0}' not found".format(refseq))
    if refproj is None:
        refproj = ""

    dproj = qsfg.DProjectNew()
    dprojref = qsfg.DProjectNew()

    dproj["name"] = os.path.splitext(output)[0]
    dproj["fName"] = output
    dproj["fNameRX"] = rxfile
    dproj["fNameRX"] = rxfile
    dproj["fNameBG"] = bgfile
    dproj["fNameSeq"] = refseq
    dproj["fNameRef"] = refproj

    if refseq is not None and refseq != "":
        dproj["RNA"] = qsff.readBaseFile(dproj["fNameSeq"])

    if refproj is not None and refproj != "":
        dprojref = getProjRefData(dproj["fNameRef"])
        dproj["RNA"] = dprojref["RNA"]
        dproj["isRef"] = True

    dataRX, dproj["Satd"]["RX"], dyesRX = qsff.readShapeData(dproj["fNameRX"])
    dataBG, dproj["Satd"]["BG"], dyesBG = qsff.readShapeData(dproj["fNameBG"])

    dproj["chIndex"] = deepcopy(channels)
    for name, idx in channels.items():
        if name.startswith("RX"):
            dproj["dData"][name] = dataRX[:, channels[name]]
            dproj["dyeN"][name] = dyesRX[channels[name]]
        else:
            dproj["dData"][name] = dataBG[:, channels[name]]
            dproj["dyeN"][name] = dyesBG[channels[name]]

    dproj["ddNTP1"] = ddNTP
    if ddNTP == "ddC":
        dproj["nuc1"] = "G"
    elif ddNTP == "ddA":
        dproj["nuc1"] = "U"
    elif ddNTP == "ddG":
        dproj["nuc1"] = "C"
    elif ddNTP == "ddT":
        dproj["nuc1"] = "A"
    else:
        raise "Invalid ddNTP"

    dproj["scriptList"] = ["New Project"]

    return dproj, dprojref


def createQuShapeFile(
    rxfile: str,
    bgfile: str,
    output: str = "output.qushape",
    refseq: str = None,
    refproj: str = None,
    channels={"RX": 0, "RXS1": 2, "BG": 0, "BGS1": 2},
    ddNTP: str = "ddC",
    path: str = None,
):

    if refseq is None and refproj is None:
        print(
            "You must provide a reference project and/or a reference sequence"
            "to generate a valid qushape project"
        )
        raise SystemExit(-1)
    dvar = qsfg.DVar(qsfg.chKeysRS)
    dproj, dprojref = createDProject(
        rxfile=rxfile,
        bgfile=bgfile,
        output=output,
        refseq=refseq,
        refproj=refproj,
        channels=channels,
        ddNTP=ddNTP,
        path=path,
    )

    if os.path.exists(output):
        os.remove(output)
    db = bsddb.hashopen(output)
    dBase = shelve.BsdDbShelf(db, protocol=2)
    dBase["dProject"] = deepcopy(dproj)
    dBase["intervalData"] = [deepcopy(dproj)]
    dBase["dProjRef"] = deepcopy(dprojref)
    dBase["dVar"] = deepcopy(dvar)
    dBase.close()
    db.close()


if __name__ == "__main__":
    fire.Fire(createQuShapeFile)
