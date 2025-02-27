import numpy as np
from datetime import datetime
import sys
import re
import glob
import os
import json
from pathlib import Path
import pandas as pd
import pickle
#this script extracts effective data from pkl files

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()

N=int(sys.argv[1])
dataRoot=f"../dataAll/N{N}/"
obs_U_dipole="U_dipole"
#search directory
TVals=[]
TFileNames=[]
TStrings=[]
for TFile in glob.glob(dataRoot+"/T*"):
    # print(TFile)
    matchT=re.search(r"T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",TFile)
    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))
        TStrings.append("T"+matchT.group(1))

#sort T values
sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]
sortedTStrings=[TStrings[ind] for ind in sortedInds]


def parseSummary(oneTFolder,obs_name):

    startingFileInd=-1

    lag=-1
    sweep_to_write=-1
    smrFile=oneTFolder+"/summary_"+obs_name+".txt"
    summaryFileExists=os.path.isfile(smrFile)
    if summaryFileExists==False:
        return startingFileInd,-1

    with open(smrFile,"r") as fptr:
        lines=fptr.readlines()

    for oneLine in lines:
        #match startingFileInd
        matchStartingFileInd=re.search(r"startingFileInd=(\d+)",oneLine)
        if matchStartingFileInd:
            startingFileInd=int(matchStartingFileInd.group(1))

        #match lag
        matchLag=re.search(r"lag=(\d+)",oneLine)
        if matchLag:
            lag=int(matchLag.group(1))
        #match sweep_to_write
        match_sweep_to_write=re.search(r"sweep_to_write=(\d+)",oneLine)

        if match_sweep_to_write:
            sweep_to_write=int(match_sweep_to_write.group(1))
    return startingFileInd,lag,sweep_to_write

def sort_data_files_by_flushEnd(oneTFolder,obs_name,varName):
    dataFolderName=oneTFolder+"/"+obs_name+"_dataFiles/"+varName+"/"
    dataFilesAll=[]
    flushEndAll=[]
    for oneDataFile in glob.glob(dataFolderName+"/flushEnd*.pkl"):
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"flushEnd(\d+)",oneDataFile)
        if matchEnd:
            flushEndAll.append(int(matchEnd.group(1)))

    endInds=np.argsort(flushEndAll)
    sortedDataFiles=[dataFilesAll[i] for i in endInds]

    return sortedDataFiles

def U_extract_ForOneT(oneTFolder,oneTStr,startingFileInd,lag):
    TRoot=oneTFolder
    sortedUDataFilesToRead=sort_data_files_by_flushEnd(TRoot,obs_U_dipole,"U")

    startingUFileName=sortedUDataFilesToRead[startingFileInd]
    with open(startingUFileName,"rb") as fptr:
        inUStart=pickle.load(fptr)

    UVec=inUStart
    for pkl_file in sortedUDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            in_UArr=pickle.load(fptr)
            UVec=np.append(UVec,in_UArr)
    UVecSelected=UVec[::lag]

    return UVecSelected


def one_dipole_component_extract_ForOneT(oneTFolder,oneTStr,startingFileInd,lag,component_name,sweep_to_write):
    TRoot=oneTFolder
    sorted_one_component_DataFilesToRead=sort_data_files_by_flushEnd(TRoot,obs_U_dipole,component_name)

    one_component_StaringFileName=sorted_one_component_DataFilesToRead[startingFileInd]

    with open(one_component_StaringFileName,"rb") as fptr:
        one_component_inArrStart=np.array(pickle.load(fptr))
    one_component_Arr=one_component_inArrStart.reshape((sweep_to_write,-1))

    #read the rest of one_v pkl files
    for pkl_file in sorted_one_component_DataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            one_component_inArr=np.array(pickle.load(fptr))
            one_component_inArr=one_component_inArr.reshape((sweep_to_write,-1))
            one_component_Arr=np.concatenate((one_component_Arr,one_component_inArr),axis=0)


    one_component_ArrSelected=one_component_Arr[::lag,:]

    return one_component_ArrSelected


def save_U_data(UVecSelected,oneTStr):
    outCsvDataRoot=dataRoot+"/csvOutAll/"
    outCsvFolder=outCsvDataRoot+"/"+oneTStr+"/"
    Path(outCsvFolder).mkdir(exist_ok=True,parents=True)

    outFileName="U.csv"
    outCsvFile=outCsvFolder+outFileName
    df=pd.DataFrame(UVecSelected)

    # Save to CSV
    df.to_csv(outCsvFile, index=False, header=False)


def save_oneComponent_dipole_data(one_component_ArrSelected,oneTStr,component_name):
    outCsvDataRoot=dataRoot+"/csvOutAll/"
    outCsvFolder=outCsvDataRoot+"/"+oneTStr+"/"
    Path(outCsvFolder).mkdir(exist_ok=True,parents=True)
    outFileName=f"{component_name}.csv"

    outCsvFile=outCsvFolder+outFileName
    df=pd.DataFrame(one_component_ArrSelected)
    # Save to CSV
    df.to_csv(outCsvFile, index=False, header=False)


for k in range(0,len(sortedTFiles)):
    tStart=datetime.now()
    oneTFolder=sortedTFiles[k]
    oneTStr=sortedTStrings[k]
    startingfileIndTmp,lagTmp,sweep_to_writeTmp=parseSummary(oneTFolder,obs_U_dipole)
    if startingfileIndTmp<0:
        print("summary file does not exist for "+oneTStr+" "+obs_U_dipole)
        continue

    UVecSelected=U_extract_ForOneT(oneTFolder,oneTStr,startingfileIndTmp,lagTmp)

    save_U_data(UVecSelected,oneTStr)
    component_Px="Px"
    component_Py="Py"
    component_Qx="Qx"
    component_Qy="Qy"

    Px_ArrSelected=one_dipole_component_extract_ForOneT(oneTFolder,oneTStr,startingfileIndTmp,lagTmp,component_Px,sweep_to_writeTmp)
    Py_ArrSelected=one_dipole_component_extract_ForOneT(oneTFolder,oneTStr,startingfileIndTmp,lagTmp,component_Py,sweep_to_writeTmp)

    Qx_ArrSelected=one_dipole_component_extract_ForOneT(oneTFolder,oneTStr,startingfileIndTmp,lagTmp,component_Qx,sweep_to_writeTmp)

    Qy_ArrSelected=one_dipole_component_extract_ForOneT(oneTFolder,oneTStr,startingfileIndTmp,lagTmp,component_Qy,sweep_to_writeTmp)

    save_oneComponent_dipole_data(Px_ArrSelected,oneTStr,component_Px)

    save_oneComponent_dipole_data(Py_ArrSelected,oneTStr,component_Py)

    save_oneComponent_dipole_data(Qx_ArrSelected,oneTStr,component_Qx)

    save_oneComponent_dipole_data(Qy_ArrSelected,oneTStr,component_Qy)

    tEnd=datetime.now()
    print("processed T="+str(sortedTVals[k])+": ",tEnd-tStart)