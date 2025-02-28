from pathlib import Path
from decimal import Decimal, getcontext

import numpy as np
import pandas as pd

#this script creates slurm bash files for exec_checking.py
def format_using_decimal(value, precision=15):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)


outPath="./bashFiles_checking/"
Path(outPath).mkdir(exist_ok=True,parents=True)
N=5 #unit cell number
TVals=[0.1,1,2,3,4,5,6,7,8,9,10]

TStrAll=[]
# print(TDirsAll)
for k in range(0,len(TVals)):
    T=TVals[k]
    # print(T)

    TStr=str(T)#format_using_decimal(T)
    TStrAll.append(TStr)


def contents_to_bash(k):
    TStr=TStrAll[k]
    contents=[
        "#!/bin/bash\n",
        "#SBATCH -n 1\n",
        "#SBATCH -N 1\n",
        "#SBATCH -t 0-60:00\n",
        "#SBATCH -p CLUSTER\n"
        "#SBATCH --mem=8GB\n",
        f"#SBATCH -o outmcT{TStr}.out\n",
        f"#SBATCH -e outmcT{TStr}.out\n",
        "cd /home/cywanag/data/hpc/cywanag/liuxi/Document/cppCode/C3_fer\n",
        f"python3 -u exec_checking.py {TStr} {N} 15\n"
        ]

    outBashName=outPath+f"/run_mcT{TStr}_checking.sh"
    with open(outBashName,"w+") as fptr:
        fptr.writelines(contents)

for k in range(0,len(TStrAll)):
    contents_to_bash(k)