## File partially from QuShape https://github.com/Weeks-UNC/QuShape
## Licenced under GPL3
## Copyright (C) 2012  Fethullah Karabiber, Oleg V. Favorov, Kevin M. Weeks


import numpy as np
#chKeysCC=['RX','RXS1','RXS2','BG','BGS1','BGS2']
chKeysRS=['RX','BG','RXS1','BGS1','RXS2','BGS2']
#chKeysRX=['RX','RXS1','RXS2']
#chKeysBG=['BG','BGS1','BGS2']

def DData(isSeq2=False):
    if isSeq2:
        keys=['RX','BG','RXS1','BGS1','RXS2','BGS2']
    else:
        keys=['RX','BG','RXS1','BGS1']
    dData={}
    for key in keys:
        dData[key]=np.array([],dtype='f4')
    return dData

def DProjectNew():
    dProjectNew={}
## NEW PROJECT RESULTS
    dProjectNew['dData']=DData() #DData(dChKeys['RS'])# {'RX':np.array([0]),'BG':np.array([0]),'RXS':np.array([0]),'BGS':np.array([0])}
    dProjectNew['dir']= ''# '' # Working directory "/Users/fethullah/example4"
    dProjectNew['name']='' #''  # Project Name
    dProjectNew['fName']='' # Project file name
    dProjectNew['fNameRX']='' # Project file name
    dProjectNew['fNameBG']='' # Project file name
    dProjectNew['fNameSeq']='' # Project file name
    dProjectNew['chIndex']={'RX':0, 'RXS1':1,'RXS2':2,'BG':0, 'BGS1':1,'BGS2':2 }
    dProjectNew['isRef']=False
    dProjectNew['fNameRef']='' # Project file name
    dProjectNew['Satd']={'RX':[],'BG':[]}
    dProjectNew['dyeN']=DData() #dict(zip(chKeys,['6-FAM','VIC','NED','6-FAM','VIC','NED']))
    dProjectNew['isSeq2']=False
    dProjectNew['RNA']=''  # Official RNA sequence
    dProjectNew['ddNTP1']='ddC' # ddNTP type of the sequences, RX BG
    dProjectNew['ddNTP2']='ddT' # ddNTP type of the sequences, RX BG
    dProjectNew['nuc1']='G'  # Nucleotide type of sequences
    dProjectNew['nuc2']='A'  # Nucleotide type of sequence
    dProjectNew['chKeyRS']=['RX', 'BG', 'RXS1', 'BGS1']
    dProjectNew['chKeyCC']=['RX', 'RXS1','BG', 'BGS1']
    dProjectNew['chKeyRX']=['RX', 'RXS1']
    dProjectNew['chKeyBG']=['BG', 'BGS1']
    dProjectNew['scriptList']=[]
    
## SEQUENCE ALIGNMENT RESULTS
    dProjectNew['seqRNA']=''
    dProjectNew['seqX0']=np.array([],dtype='i4') # X value of nuc, above threshold 
    dProjectNew['seq0']='' #Detected sequence
    dProjectNew['seqX']=np.array([],dtype='i4') # X value of completed seq
    dProjectNew['usedCh1']='BGS1' # used channel RXS1 or BGS1
    dProjectNew['usedCh2']='BGS2' # used channel RXS1 or BGS1
    dProjectNew['start']=0  # start point in RNA, roi 
    dProjectNew['end']=0# end point in RNA, roi 
    dProjectNew['seqNum']=np.array([],int)# end point in RNA, roi 
    dProjectNew['scrNuc']=np.array([],dtype='f4')
    
## REACTIVITY RESULTS 
    dProjectNew['dPeakRX']=DPeakList()
    dProjectNew['dPeakBG']=DPeakList()
 #   dProjectNew['dPeakBG']=DPeakList()
    
    dProjectNew['scaleFactor']=np.array([1])
    dProjectNew['areaDiff']=np.array([],dtype='f4')
    dProjectNew['normDiff']=np.array([],dtype='f4')
    
  
    return dProjectNew

def DPeakList():
    dPeakList={}
    dPeakList['NPeak']=0
    dPeakList['pos']=np.array([],dtype='i4')
    dPeakList['amp']=np.array([],dtype='f4')
    dPeakList['wid']=np.array([],dtype='f4')
    dPeakList['area']=np.array([],dtype='f4')
    dPeakList['averW']=np.array([],dtype='f4')
    dPeakList['minW']=np.array([],dtype='f4')
    dPeakList['maxW']=np.array([],dtype='f4')
    
    return dPeakList

def DVar(chKeyRS):
    colors=['#FF0000','#0000FF','#008000',
            '#FF00FF','#FFA500','#00FFFF']
    dVar={}
    dVar['lineVisible']={}
    dVar['lineColor']={}
    dVar['lineStyle']={}
    dVar['lineMarker']={}
    dVar['lineWidth']={}
    
    for i in range(len(chKeyRS)):
        key=chKeyRS[i]
        dVar['lineColor'][key]=colors[i]  
        dVar['lineVisible'][key]=True
        dVar['lineStyle'][key]='-'
        dVar['lineMarker'][key]=''
        dVar['lineWidth'][key]=1
    
    dVar['widthP']=100
    dVar['heightP']=100
    dVar['zoomP']=100
    dVar['left']=0.01
    dVar['right']=0.99
    dVar['top']=0.99
    dVar['bottom']=0.05
    
    dVar['maxLength']=0

#
    dVar['isDoneSeqAlign']=False
    dVar['isDoneSeqAlignRef']=False
    dVar['isDoneReactivity']=False
    
## CONTROL FLAGS
    dVar['flag']={}
    dVar['flag']['isDrawStad']=False
    dVar['flag']['isSeqAlign']=False
    dVar['flag']['isDrawLine']=False
    dVar['flag']['isDrawGauss']=False
    
    dVar['flag']['isSatd']=False
    dVar['flag']['isScale']=False
    dVar['flag']['isPeakMatchModify']=False
    dVar['flag']['isPeakLinkRefModify']=False
    dVar['flag']['isDrawRef']=False
    
    
    return dVar
    
def globalVars():
    dGlobalVars={}
    dGlobalVars['dye']={'5-FAM':517,'6-FAM':522,'TET':538,'HEX':553,'JOE':554,'VIC':555,'NED':575,'TAMRA':583,'PET':595,'ROX':607,'LIZ':655}
    return dGlobalVars
       
