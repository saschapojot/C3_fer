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


#this script concatenates pkl files and plot U
def format_using_decimal(value, precision=15):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)
T=4
N=5
startingFileInd=40
TStr=format_using_decimal(T)
dataDirRoot=f"./dataAll/N{N}/T{TStr}/U_dipole_dataFiles/"


var_name="/U/"

dataIn_dir=dataDirRoot+var_name


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


sortedDataFiles=sort_data_files_by_flushEnd(dataIn_dir)

with open(sortedDataFiles[startingFileInd],"rb")as fptr:
    inArrStart=pickle.load(fptr)
print(f"len(sortedDataFiles)={len(sortedDataFiles)}")
U_arr=inArrStart
for pkl_file in sortedDataFiles[startingFileInd+1:]:
    with open(pkl_file,"rb") as fptr:
        inArr=pickle.load(fptr)
    U_arr=np.append(U_arr,inArr)

out_pic_root=f"./dataAll/N{N}/"
plt.figure()
plt.plot(range(0,len(U_arr)),U_arr,color="blue")
plt.ylabel("$U$")
plt.title(f"T={TStr}, startingFileInd={startingFileInd}")
plt.savefig(out_pic_root+f"tmpU_N{N}_T{TStr}.png")
plt.close()
