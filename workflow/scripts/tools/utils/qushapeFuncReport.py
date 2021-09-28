## File partially from QuShape https://github.com/Weeks-UNC/QuShape
## Licenced under GPL3
## Copyright (C) 2012  Fethullah Karabiber, Oleg V. Favorov, Kevin M. Weeks

import numpy as np


reportKeys=['seqNum','seqRNA','posSeq','posRX','areaRX','posBG','areaBG','areaDiff','normDiff']
def DReport():
    dReport={}
    dReport['seqNum']=np.array([],dtype='i4')
    dReport['seqRNA']=''
    dReport['posSeq']=np.array([],dtype='i4')
    dReport['posRX']=np.array([],dtype='i4')
    dReport['areaRX']=np.array([],dtype='i4')
    dReport['posBG']=np.array([],dtype='i4')
    dReport['areaBG']=np.array([],dtype='i4')
    dReport['areaDiff']=np.array([],dtype='i4')
    dReport['normDiff']=np.array([],dtype='i4')
    
    return dReport
 
def createDReport(dProject):
    dReport=DReport()
    dReport['seqNum']=np.array(dProject['seqNum'][1:],int)
    dReport['seqRNA']=dProject['seqRNA'][1:]
    dReport['posSeq']=np.array(dProject['seqX'][1:],int)
    dReport['posRX']=np.array(dProject['dPeakRX']['pos'][:-1],int)
    dReport['areaRX']=np.round(dProject['dPeakRX']['area'][:-1],decimals=2)
    dReport['posBG']=np.array(dProject['dPeakBG']['pos'][:-1],int)
    dReport['areaBG']=np.round(dProject['dPeakBG']['area'][:-1],decimals=2)
    dReport['areaDiff']=np.round(dProject['areaDiff'][:-1],decimals=2)
    dReport['normDiff']=np.round(dProject['normDiff'][:-1],decimals=2)
    return dReport

def writeReportFile(dReport,fName):
    myfile=open(fName,'w')    
    for key in reportKeys:
        myfile.write(str(key)+'\t')
    myfile.write('\n')
    for i in range(len(dReport['seqRNA'])):
        for key in reportKeys:
            myfile.write(str(dReport[key][i])+'\t')
        myfile.write('\n')
 
