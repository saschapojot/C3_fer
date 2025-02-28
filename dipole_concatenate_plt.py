import pickle

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from multiprocessing import Pool
import pandas as pd
import statsmodels.api as sm
import sys
import re
import warnings
from decimal import Decimal, getcontext

from scipy.stats import ks_2samp
import glob

import os
import json
import pickle

#this script concatenates pkl files and plot dipole magnitude, angle

def format_using_decimal(value, precision=15):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

T=8
N=5
startingFileInd=60
TStr=format_using_decimal(T)
dataDirRoot=f"./dataAll/N{N}/T{TStr}/U_dipole_dataFiles/"
eps_for_auto_corr=1e-1
sweep_to_write=500

var_name_Px="Px"
var_name_Py="Py"
var_name_Qx="Qx"
var_name_Qy="Qy"

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
def load_one_variable(var_name):
    dataIn_dir=dataDirRoot+"/"+var_name+"/"
    sortedDataFiles=sort_data_files_by_flushEnd(dataIn_dir)

    with open(sortedDataFiles[startingFileInd],"rb")as fptr:
        inArrStart=np.array(pickle.load(fptr))

    arr=inArrStart.reshape((sweep_to_write,-1))

    for pkl_file in sortedDataFiles[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            inArr=np.array(pickle.load(fptr))
        inArr=inArr.reshape((sweep_to_write,-1))
        arr=np.concatenate((arr,inArr),axis=0)

    return arr



Px_arr=load_one_variable(var_name_Px)
Py_arr=load_one_variable(var_name_Py)
Qx_arr=load_one_variable(var_name_Qx)
Qy_arr=load_one_variable(var_name_Qy)

dipole_x=Px_arr+Qx_arr

dipole_y=Py_arr+Qy_arr

dipole_x_avg=np.mean(dipole_x,axis=1)
dipole_y_avg=np.mean(dipole_y,axis=1)

dipole_abs_avg=np.sqrt(dipole_x_avg**2+dipole_y_avg**2)

sameTmp_dipole,lagTmp_dipole=auto_corrForOneVec(dipole_abs_avg)

print(f"dipole magitude lag={lagTmp_dipole}")

out_pic_root=f"./dataAll/N{N}/"
plt.figure()
plt.plot(range(0,len(dipole_abs_avg)),dipole_abs_avg,color="green")
plt.ylabel("|dipole|")
plt.title(f"T={TStr}, startingFileInd={startingFileInd}")
plt.savefig(out_pic_root+f"tmp_dipole_abs_N{N}_T{TStr}.png")
plt.close()