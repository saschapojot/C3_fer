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

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()

N=int(sys.argv[1])
TStr=sys.argv[2]

dataRoot=f"../dataAll/N{N}/T{TStr}/"
obs_U_dipole="polarization"

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


def sort_data_files_by_flushEnd(oneTFolder,varName):
    dataFolderName=oneTFolder+"/U_dipole_dataFiles/"+varName+"/"
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

def one_dipole_component_extract_ForOneT(oneTFolder,startingFileInd,lag,component_name,sweep_to_write):
    TRoot=oneTFolder
    sorted_one_component_DataFilesToRead=sort_data_files_by_flushEnd(TRoot,component_name)

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



def save_oneComponent_dipole_data(one_component_ArrSelected,oneTStr,component_name):
    outCsvDataRoot=dataRoot+"/csvOutAll/"
    outCsvFolder=outCsvDataRoot+"/"+oneTStr+"/"

    Path(outCsvFolder).mkdir(exist_ok=True,parents=True)
    outFileName=f"{component_name}.csv"

    outCsvFile=outCsvFolder+outFileName
    df=pd.DataFrame(one_component_ArrSelected)
    # Save to CSV
    df.to_csv(outCsvFile, index=False, header=False)



startingFileInd,lag,sweep_to_write=parseSummary(dataRoot,obs_U_dipole)
if startingFileInd<0:
    print("summary file does not exist for "+TStr+" "+obs_U_dipole)
    exit(1)

component_Px="Px"
component_Py="Py"
component_Qx="Qx"
component_Qy="Qy"

Px_ArrSelected=one_dipole_component_extract_ForOneT(dataRoot,startingFileInd,lag,component_Px,sweep_to_write)

Py_ArrSelected=one_dipole_component_extract_ForOneT(dataRoot,startingFileInd,lag,component_Py,sweep_to_write)

Qx_ArrSelected=one_dipole_component_extract_ForOneT(dataRoot,startingFileInd,lag,component_Qx,sweep_to_write)


Qy_ArrSelected=one_dipole_component_extract_ForOneT(dataRoot,startingFileInd,lag,component_Qy,sweep_to_write)


save_oneComponent_dipole_data(Px_ArrSelected,TStr,component_Px)

save_oneComponent_dipole_data(Py_ArrSelected,TStr,component_Py)


save_oneComponent_dipole_data(Qx_ArrSelected,TStr,component_Qx)

save_oneComponent_dipole_data(Qy_ArrSelected,TStr,component_Qy)