import numpy as np
from datetime import datetime
import statsmodels.api as sm
import sys
import re
import warnings


from scipy.stats import ks_2samp
import glob

import os
import json
import pickle



#This script checks if U, Px,Py,Qx,Qy values reach equilibrium and writes summary file of dist
#This file checks pkl files

argErrCode=2
sameErrCode=3
missingErrCode=4
eps_for_auto_corr=5e-2
if (len(sys.argv)!=4):
    print("wrong number of arguments")
    exit(argErrCode)

jsonFromSummaryLast=json.loads(sys.argv[1])
jsonDataFromConf=json.loads(sys.argv[2])
lastFileNum=int(sys.argv[3])
TDirRoot=jsonFromSummaryLast["TDirRoot"]
U_dipole_dataDir=jsonFromSummaryLast["U_dipole_dataDir"]
effective_data_num_required=int(jsonDataFromConf["effective_data_num_required"])

N=int(jsonDataFromConf["N"])
sweep_to_write=int(jsonDataFromConf["sweep_to_write"])
summary_U_dipoleFile=TDirRoot+"/summary_U_dipole.txt"
# lastFileNum=7
def sort_data_files_by_flushEnd(oneDir):
    dataFilesAll=[]
    flushEndAll=[]
    for oneDataFile in glob.glob(oneDir+"/flushEnd*.pkl"):
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"flushEnd(\d+)",oneDataFile)
        if matchEnd:
            indTmp=int(matchEnd.group(1))
            flushEndAll.append(indTmp)

    endInds=np.argsort(flushEndAll)
    sortedDataFiles=[dataFilesAll[i] for i in endInds]
    return sortedDataFiles


def parseSummaryU_dipole():
    startingFileInd=-1
    summaryFileExists=os.path.isfile(summary_U_dipoleFile)
    if summaryFileExists==False:
        return startingFileInd
    with open(summary_U_dipoleFile,"r") as fptr:
        lines=fptr.readlines()

    for oneLine in lines:
        #match startingFileInd
        matchStartingFileInd=re.search(r"startingFileInd=(\d+)",oneLine)
        if matchStartingFileInd:
            startingFileInd=int(matchStartingFileInd.group(1))

    return startingFileInd


def auto_corrForOneVec(vec):
    """

    :param colVec: a vector of data
    :return:
    """
    same=False
    NLags=int(len(vec)*3/4)
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVec=sm.tsa.acf(vec,nlags=NLags)
    except Warning as w:
        same=True

    acfOfVecAbs=np.abs(acfOfVec)
    minAutc=np.min(acfOfVecAbs)
    lagVal=-1
    if minAutc<=eps_for_auto_corr:
        lagVal=np.where(acfOfVecAbs<=eps_for_auto_corr)[0][0]

    return same,lagVal

def ksTestOneColumn(vec,lag):
    """

    :param vec: a vector of data
    :param lag: auto-correlation length
    :return:
    """
    vecSelected=vec[::lag]
    lengthTmp=len(vecSelected)
    if lengthTmp%2==1:
        lengthTmp-=1
    lenPart=int(lengthTmp/2)
    vecToCompute=vecSelected[-lengthTmp:]

    #ks test
    selectedVecPart0=vecToCompute[:lenPart]
    selectedVecPart1=vecToCompute[lenPart:]
    result=ks_2samp(selectedVecPart0,selectedVecPart1)

    return result.pvalue,result.statistic, lenPart*2

def checkUDataFilesForOneT(UData_dir):
    U_sortedDataFilesToRead=sort_data_files_by_flushEnd(UData_dir)
    if len(U_sortedDataFilesToRead)==0:
        print("no data for U.")
        exit(0)

    startingFileInd=parseSummaryU_dipole()
    if startingFileInd<0:
        #we guess that the equilibrium starts at this file
        startingFileInd=len(U_sortedDataFilesToRead)-lastFileNum
    startingFileName=U_sortedDataFilesToRead[startingFileInd]
    with open(startingFileName,"rb") as fptr:
        inArrStart=pickle.load(fptr)

    U_arr=inArrStart
    #read the rest of the pkl files
    for pkl_file in U_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            inArr=pickle.load(fptr)
        U_arr=np.append(U_arr,inArr)

    sameUTmp,lagUTmp=auto_corrForOneVec(U_arr)

    if sameUTmp==True or lagUTmp==-1:
        return [sameUTmp,lagUTmp,-1,-1,-1,-1]

    pUTmp,statUTmp,lengthUTmp=ksTestOneColumn(U_arr,lagUTmp)
    numDataPoints=lengthUTmp

    return [sameUTmp,lagUTmp,pUTmp,statUTmp,numDataPoints,startingFileInd]



def check_DipoleDataFilesForOneT(Px_dir,Py_dir,Qx_dir,Qy_dir):

    Px_sortedDataFilesToRead=sort_data_files_by_flushEnd(Px_dir)
    # print(f"Px_sortedDataFilesToRead={Px_sortedDataFilesToRead}")
    Py_sortedDataFilesToRead=sort_data_files_by_flushEnd(Py_dir)

    Qx_sortedDataFilesToRead=sort_data_files_by_flushEnd(Qx_dir)

    Qy_sortedDataFilesToRead=sort_data_files_by_flushEnd(Qy_dir)

    len_Px=len(Px_sortedDataFilesToRead)

    len_Py=len(Py_sortedDataFilesToRead)

    len_Qx=len(Qx_sortedDataFilesToRead)

    len_Qy=len(Qy_sortedDataFilesToRead)

    diff2=(len_Px-len_Py)**2+(len_Px-len_Qx)**2+(len_Px-len_Qy)**2
    if diff2>0:
        print(f"diff2={diff2}, data missing.")
        exit(missingErrCode)
    startingFileInd=parseSummaryU_dipole()
    if startingFileInd<0:
        #we guess that the equilibrium starts at this file
        startingFileInd=len_Px-lastFileNum

    Px_startingFileName=Px_sortedDataFilesToRead[startingFileInd]

    Py_startingFileName=Py_sortedDataFilesToRead[startingFileInd]

    Qx_startingFileName=Qx_sortedDataFilesToRead[startingFileInd]

    Qy_startingFileName=Qy_sortedDataFilesToRead[startingFileInd]
    # print(f"Px_startingFileName={Px_startingFileName}")
    # print(f"Py_startingFileName={Py_startingFileName}")
    # print(f"Qx_startingFileName={Qx_startingFileName}")
    # print(f"Qy_startingFileName={Qy_startingFileName}")

    #read Px
    with open(Px_startingFileName,"rb") as fptr:
        Px_inArrStart=np.array(pickle.load(fptr))
    Px_arr=Px_inArrStart.reshape((sweep_to_write,-1))

    #read the rest of Px pkl files
    for pkl_file in Px_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            Px_inArr=np.array(pickle.load(fptr))
        Px_inArr=Px_inArr.reshape((sweep_to_write,-1))
        Px_arr=np.concatenate((Px_arr,Px_inArr),axis=0)

    #read Py
    with open(Py_startingFileName,"rb") as fptr:
        Py_inArrStart=np.array(pickle.load(fptr))

    Py_arr=Py_inArrStart.reshape((sweep_to_write,-1))
    #read the rest of Py pkl files
    for pkl_file in Py_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            Py_inArr=np.array(pickle.load(fptr))
        Py_inArr=Py_inArr.reshape((sweep_to_write,-1))
        Py_arr=np.concatenate((Py_arr,Py_inArr),axis=0)

    #read Qx
    with open(Qx_startingFileName,"rb") as fptr:
        Qx_inArrStart=np.array(pickle.load(fptr))
    Qx_arr=Qx_inArrStart.reshape((sweep_to_write,-1))
    #read the rest of Qx pkl files
    for pkl_file in Qx_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            Qx_inArr=np.array(pickle.load(fptr))
        Qx_inArr=Qx_inArr.reshape((sweep_to_write,-1))
        Qx_arr=np.concatenate((Qx_arr,Qx_inArr),axis=0)

    #read Qy
    with open(Qy_startingFileName,"rb") as fptr:
        Qy_inArrStart=np.array(pickle.load(fptr))
    Qy_arr=Qy_inArrStart.reshape((sweep_to_write,-1))
    #read the rest of Qy pkl files
    for pkl_file in Qy_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            Qy_inArr=np.array(pickle.load(fptr))
        Qy_inArr=Qy_inArr.reshape((sweep_to_write,-1))
        Qy_arr=np.concatenate((Qy_arr,Qy_inArr),axis=0)

    Px_avg=np.mean(Px_arr,axis=1)

    Py_avg=np.mean(Py_arr,axis=1)

    Qx_avg=np.mean(Qx_arr,axis=1)

    Qy_avg=np.mean(Qy_arr,axis=1)
    # print(Qy_avg[-40:])



    sameTmp_Px,lagTmp_Px=auto_corrForOneVec(Px_avg)

    sameTmp_Py,lagTmp_Py=auto_corrForOneVec(Py_avg)

    sameTmp_Qx,lagTmp_Qx=auto_corrForOneVec(Qx_avg)

    sameTmp_Qy,lagTmp_Qy=auto_corrForOneVec(Qy_avg)

    same_vec_PxPyQxQy=[sameTmp_Px,sameTmp_Py,sameTmp_Qx,sameTmp_Qy]
    lag_vec_PxPyQxQy=[lagTmp_Px,lagTmp_Py,lagTmp_Qx,lagTmp_Qy]
    # print(f"lag_vec_PxPyQxQy={lag_vec_PxPyQxQy}")
    if any(same_vec_PxPyQxQy) or -1 in lag_vec_PxPyQxQy:
        return [-2],[-1],-1,[]



    lagMax=np.max(lag_vec_PxPyQxQy)

    pTmp_Px,statTmp_Px,lengthTmp_Px=ksTestOneColumn(Px_avg,lagMax)

    pTmp_Py,statTmp_Py,lengthTmp_Py=ksTestOneColumn(Py_avg,lagMax)

    pTmp_Qx,statTmp_Qx,lengthTmp_Qx=ksTestOneColumn(Qx_avg,lagMax)

    pTmp_Qy,statTmp_Qy,lengthTmp_Qy=ksTestOneColumn(Qy_avg,lagMax)

    pVec=[pTmp_Px,pTmp_Py,pTmp_Qx,pTmp_Qy]
    statVec=[statTmp_Px,statTmp_Py,statTmp_Qx,statTmp_Qy]

    numDataPoints=lengthTmp_Qy
    return pVec,statVec,numDataPoints,lag_vec_PxPyQxQy






UDataDir=U_dipole_dataDir+"/U/"
sameVec=[]
lagVec=[]
pVec=[]
statVec=[]
numDataVec=[]
print("checking U")
sameUTmp,lagUTmp,pUTmp,statUTmp,numDataPointsU,startingFileInd=checkUDataFilesForOneT(UDataDir)
sameVec.append(sameUTmp)
lagVec.append(lagUTmp)
pVec.append(pUTmp)
statVec.append(statUTmp)
numDataVec.append(numDataPointsU)
print("lagU="+str(lagUTmp))


Px_dir=U_dipole_dataDir+"/Px/"
Py_dir=U_dipole_dataDir+"/Py/"
Qx_dir=U_dipole_dataDir+"/Qx/"
Qy_dir=U_dipole_dataDir+"/Qy/"

pVec_vec_PxPyQxQy,statVec_vec_PxPyQxQy,numDataPoints_PxPyQxQy,lag_vec_PxPyQxQy=check_DipoleDataFilesForOneT(Px_dir,Py_dir,Qx_dir,Qy_dir)

pVec+=pVec_vec_PxPyQxQy
statVec+=statVec_vec_PxPyQxQy

lagVecAll=[lagUTmp]+lag_vec_PxPyQxQy
lagMax=np.max(lagVecAll)
numDataPoints=np.min([numDataPointsU,numDataPoints_PxPyQxQy])
print("lagMax="+str(lagMax))
print("numDataPoints="+str(numDataPoints))
print(f"pVec={pVec}")
print(f"statVec={statVec}")
############################################

statThreshhold=0.1
if pVec[0]==-2:
    with open(summary_U_dipoleFile,"w+") as fptr:
        msg="error: same\n"
        fptr.writelines(msg)
        exit(sameErrCode)
if numDataPoints<0:
    msg="high correlation"
    with open(summary_U_dipoleFile,"w+") as fptr:
        fptr.writelines(msg)
    exit(0)

if (np.min(pVec)>=0.01 or np.max(statVec)<=statThreshhold) and numDataPoints>=200:
    if numDataPoints>=effective_data_num_required:
        newDataPointNum=0
    else:
        newDataPointNum=effective_data_num_required-numDataPoints

    msg="equilibrium\n" \
        +"lag="+str(lagMax)+"\n" \
        +"numDataPoints="+str(numDataPoints)+"\n" \
        +"startingFileInd="+str(startingFileInd)+"\n" \
        +"newDataPointNum="+str(newDataPointNum)+"\n" \
        +"sweep_to_write="+str(sweep_to_write)+"\n"

    print(msg)
    with open(summary_U_dipoleFile,"w+") as fptr:
        fptr.writelines(msg)
    exit(0)
#continue
continueMsg="continue\n"
if not (np.min(pVec)>=0.01 or np.max(statVec)<=statThreshhold):
    continueMsg+="stat value: "+str(np.max(statVec))+"\n"
    continueMsg+="p value: "+str(np.min(pVec))+"\n"


if numDataPoints<200:
    #not enough data number

    continueMsg+="numDataPoints="+str(numDataPoints)+" too low\n"
    continueMsg+="lag="+str(lagMax)+"\n"
print(continueMsg)
with open(summary_U_dipoleFile,"w+") as fptr:
    fptr.writelines(continueMsg)
exit(0)
