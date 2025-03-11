import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd
import scipy.stats as stats

#This script loads csv data and plot P, with confidence interval


if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()

N=int(sys.argv[1])
csvDataFolderRoot=f"../dataAll/N{N}/csvOutAll/"
TVals=[]
TFileNames=[]

for TFile in glob.glob(csvDataFolderRoot+"/T*"):

    matchT=re.search(r"T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",TFile)
    # if float(matchT.group(1))<1:
    #     continue

    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))

sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]


def vec_estimator(vec):
    E_vec=np.mean(vec)
    return E_vec


def vec_jackknife(vec):
    n=len(vec)
    jackknife_samples = np.zeros(n)
    for i  in range(0,n):
        sample_vec=np.delete(vec,i)
        jackknife_samples[i]=vec_estimator(sample_vec)
    # Jackknife estimate of the statistic
    jackknife_estimate = np.mean(jackknife_samples)
    variance_estimate = (n - 1) / n * np.sum((jackknife_samples - jackknife_estimate) ** 2)

    return jackknife_estimate, variance_estimate

def vec_confidence_interval(vec,confidence_level=0.95):
    jackknife_estimate, jackknife_variance=vec_jackknife(vec)

    n=len(vec)
    alpha = 1 - confidence_level
    t_critical = stats.t.ppf(1 - alpha / 2, df=n-1)
    # Calculate the standard error
    standard_error = np.sqrt(jackknife_variance)
    # Calculate the confidence interval
    ci_lower = jackknife_estimate - t_critical * standard_error
    ci_upper = jackknife_estimate + t_critical * standard_error
    return jackknife_estimate,ci_lower, ci_upper


def generate_one_point_P_abs(oneTFile):
    matchT=re.search(r'T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',oneTFile)
    TVal=float(matchT.group(1))

    Px_path=oneTFile+"/Px.csv"
    Py_path=oneTFile+"/Py.csv"

    df_Px=np.array(pd.read_csv(Px_path,header=None))
    df_Py=np.array(pd.read_csv(Py_path,header=None))

    Px_avg_vec=np.mean(df_Px,axis=1)
    Py_avg_vec=np.mean(df_Py,axis=1)

    P_abs_vec=np.sqrt(Px_avg_vec**2+Py_avg_vec**2)
    print(f"T={TVal}, data num={len(P_abs_vec)}")
    jackknife_estimate,ci_lower, ci_upper=vec_confidence_interval(P_abs_vec)

    return [jackknife_estimate,ci_lower, ci_upper]



P_abs_valsAll=[]
interval_lowerValsAll=[]
interval_upperValsAll=[]
for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    jackknife_estimate,ci_lower, ci_upper=generate_one_point_P_abs(oneTFile)
    P_abs_valsAll.append(jackknife_estimate)
    interval_lowerValsAll.append(ci_lower)
    interval_upperValsAll.append(ci_upper)

sortedTVals=np.array(sortedTVals)
TInds=np.where(sortedTVals>0.2)
TToPlt=sortedTVals[TInds]
interval_lowerValsAll=np.array(interval_lowerValsAll)
interval_upperValsAll=np.array(interval_upperValsAll)
P_abs_valsAll=np.array(P_abs_valsAll)

P_abs_br=P_abs_valsAll-interval_lowerValsAll
fig,ax=plt.subplots()
ax.errorbar(TToPlt,P_abs_valsAll[TInds],yerr=P_abs_br[TInds],fmt='o',color="black", ecolor='r', capsize=5,label='mc')
ax.set_xscale("log")
ax.set_xlabel('$T$')
ax.set_ylabel("$|P|$")
ax.set_title("$|P|$ per unit cell, unit cell number="+str(N**2))
plt.legend(loc="best")
plt.savefig(csvDataFolderRoot+"/P_absPerUnitCell.png")
plt.close()
