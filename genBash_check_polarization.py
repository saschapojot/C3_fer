from pathlib import Path
from decimal import Decimal, getcontext
import shutil
import numpy as np
import pandas as pd
import os
#this script creates slurm bash files for check_polarization_after_one_run.py
def format_using_decimal(value, precision=4):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

outPath="./bashFiles_checking/"
if os.path.isdir(outPath):
    shutil.rmtree(outPath)

Path(outPath).mkdir(exist_ok=True,parents=True)
N=5 #unit cell number
number=int((10-0.1)/0.01)
TVals=[0.1+0.01*n for n in range(0,number+1)]
lastFileNum=20
TStrAll=[]
chunk_size = 100
chunks = [TVals[i:i + chunk_size] for i in range(0, len(TVals), chunk_size)]

def contents_to_bash(chk_ind,T_ind,chunks):

    TStr=format_using_decimal(chunks[chk_ind][T_ind])
    conf_file_name=f"./dataAll/N{N}/T{TStr}/run_T{TStr}.mc.conf"
    contents=[
        "#!/bin/bash\n",
        "#SBATCH -n 2\n",
        "#SBATCH -N 1\n",
        "#SBATCH -t 0-60:00\n",
        "#SBATCH -p hebhcnormal01\n",
        "#SBATCH --mem=4GB\n",
        f"#SBATCH -o out_polarization_{TStr}.out\n",
        f"#SBATCH -e out_polarization_{TStr}.err\n",
        "cd /public/home/hkust_jwliu_1/liuxi/Document/cppCode/C3_fer\n",
        f"python3 -u check_polarization_after_one_run.py {conf_file_name} {lastFileNum}\n"
        ]

    out_chunk=outPath+f"/chunk{chk_ind}/"
    Path(out_chunk).mkdir(exist_ok=True,parents=True)
    outBashName=out_chunk+f"/check_pol_T{TStr}.sh"
    with open(outBashName,"w+") as fptr:
        fptr.writelines(contents)


for chk_ind in range(0,len(chunks)):
    for T_ind in range(0,len(chunks[chk_ind])):
        contents_to_bash(chk_ind,T_ind,chunks)