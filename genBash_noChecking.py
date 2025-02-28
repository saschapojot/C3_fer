from pathlib import Path
from decimal import Decimal, getcontext
import os
import numpy as np
import pandas as pd

#this script creates slurm bash files for exec_noChecking.py
def format_using_decimal(value, precision=15):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)


outPath="./bashFiles_noChecking/"
os.rmdir(outPath)
Path(outPath).mkdir(exist_ok=True,parents=True)
N=5 #unit cell number
TVals=[10,11,13,15]

TStrAll=[]
for k in range(0,len(TVals)):
    T=TVals[k]
    # print(T)

    TStr=str(T)#format_using_decimal(T)
    TStrAll.append(TStr)

def contents_to_bash_noChecking(k):
    TStr=TStrAll[k]
    contents=[
        "#!/bin/bash\n",
        "#SBATCH -n 1\n",
        "#SBATCH -N 1\n",
        "#SBATCH -t 0-60:00\n",
        "#SBATCH -p CLUSTER\n"
        "#SBATCH --mem=8GB\n",
        f"#SBATCH -o outmcT{TStr}.out\n",
        f"#SBATCH -e outmcT{TStr}.err\n",
        "cd /home/cywanag/data/hpc/cywanag/liuxi/Document/cppCode/C3_fer\n",
        f"python3 -u exec_noChecking.py {TStr} {N}\n"
    ]

    outBashName=outPath+f"/run_mcT{TStr}_noChecking.sh"
    with open(outBashName,"w+") as fptr:
        fptr.writelines(contents)



for k in range(0,len(TStrAll)):
    contents_to_bash_noChecking(k)