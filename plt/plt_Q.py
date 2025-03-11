import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd
import scipy.stats as stats

#This script loads csv data and plot Q, with confidence interval

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
    avg_vec=np.mean(vec)
    return avg_vec

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


def generate_one_point_Q_abs(oneTFile):
    matchT=re.search(r'T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',oneTFile)
    TVal=float(matchT.group(1))

    Qx_path=oneTFile+"/Qx.csv"
    Qy_path=oneTFile+"/Qy.csv"

    df_Qx=np.array(pd.read_csv(Qx_path,header=None))
    df_Qy=np.array(pd.read_csv(Qy_path,header=None))

    Qx_avg_vec=np.mean(df_Qx,axis=1)
    Qy_avg_vec=np.mean(df_Qy,axis=1)

    Q_abs_vec=np.sqrt(Qx_avg_vec**2+Qy_avg_vec**2)
    print(f"T={TVal}, data num={len(Q_abs_vec)}")

    jackknife_estimate,ci_lower, ci_upper=vec_confidence_interval(Q_abs_vec)

    return [jackknife_estimate,ci_lower, ci_upper]



Q_abs_valsAll=[]
interval_lowerValsAll=[]
interval_upperValsAll=[]
for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    jackknife_estimate,ci_lower, ci_upper=generate_one_point_Q_abs(oneTFile)
    Q_abs_valsAll.append(jackknife_estimate)
    interval_lowerValsAll.append(ci_lower)
    interval_upperValsAll.append(ci_upper)


sortedTVals=np.array(sortedTVals)
TInds=np.where(sortedTVals>0.2)
TToPlt=sortedTVals[TInds]
interval_lowerValsAll=np.array(interval_lowerValsAll)
interval_upperValsAll=np.array(interval_upperValsAll)

Q_abs_valsAll=np.array(Q_abs_valsAll)
Q_abs_br=Q_abs_valsAll-interval_lowerValsAll
fig,ax=plt.subplots()
ax.errorbar(TToPlt,Q_abs_valsAll[TInds],yerr=Q_abs_br[TInds],fmt='o',color="black", ecolor='r', capsize=5,label='mc')
ax.set_xlabel('$T$')
ax.set_xscale("log")
ax.set_ylabel("$|Q|$")
ax.set_title("$|Q|$ per unit cell, unit cell number="+str(N**2))
plt.legend(loc="best")
plt.savefig(csvDataFolderRoot+"/Q_absPerUnitCell.png")
plt.close()