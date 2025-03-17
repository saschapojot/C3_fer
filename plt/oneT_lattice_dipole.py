import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd
import scipy.stats as stats
import os
from decimal import Decimal, getcontext
errFileNotExist=1
def format_using_decimal(value, precision=4):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)
#this script converts Px, Py, Qx, Qy csv files to average

#This script loads csv data and plot P+Q, with confidence interval
if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()

N=int(sys.argv[1])
TStr=sys.argv[2]
# csvDataFolderRoot=f"../dataAll/N{N}/csvOutAll/T{TStr}/"
print(f"csvDataFolderRoot={csvDataFolderRoot}")
Px_csv_file=csvDataFolderRoot+"/Px.csv"

Py_csv_file=csvDataFolderRoot+"/Py.csv"

Qx_csv_file=csvDataFolderRoot+"/Qx.csv"

Qy_csv_file=csvDataFolderRoot+"/Qy.csv"

if not os.path.exists(Px_csv_file):
    print(f"Px does not exist for {TStr}")
    exit(errFileNotExist)

if not os.path.exists(Py_csv_file):
    print(f"Py does not exist for {TStr}")
    exit(errFileNotExist)


if not os.path.exists(Qx_csv_file):
    print(f"Qx does not exist for {TStr}")
    exit(errFileNotExist)


if not os.path.exists(Qy_csv_file):
    print(f"Qy does not exist for {TStr}")
    exit(errFileNotExist)

Px_csv_arr=np.array(pd.read_csv(Px_csv_file,header=None))

Px_avg=np.mean(Px_csv_arr,axis=0)

Py_csv_arr=np.array(pd.read_csv(Py_csv_file,header=None))

Py_avg=np.mean(Py_csv_arr,axis=0)

Qx_csv_arr=np.array(pd.read_csv(Qx_csv_file,header=None))

Qx_avg=np.mean(Qx_csv_arr,axis=0)

Qy_csv_arr=np.array(pd.read_csv(Qy_csv_file,header=None))

Qy_avg=np.mean(Qy_csv_arr,axis=0)

out_dipole_file_name=csvDataFolderRoot+"/avg_dipole_combined.csv"

out_arr=np.array([
    Px_avg,Py_avg,Qx_avg,Qy_avg
])

df=pd.DataFrame(out_arr)
df.to_csv(out_dipole_file_name, header=False, index=False)




