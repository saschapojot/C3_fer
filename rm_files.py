import subprocess
import signal
import sys
import os
import re
import glob
from datetime import datetime
import numpy as np

import shutil

#this script removes some files

if len(sys.argv) != 2:
    print("wrong number of arguments")
    sys.exit(1)

N = int(sys.argv[1])
dataRoot = f"./dataAll/N{N}/csvOutAll/"

# search directory
TVals = []
TFileNames = []
TStrings = []

for TFile in glob.glob(os.path.join(dataRoot, "T*")):
    matchT = re.search(r"T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)", TFile)
    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))
        TStrings.append(matchT.group(1))


# sort T values
sortedInds = np.argsort(TVals)
sortedTVals = [TVals[ind] for ind in sortedInds]
sortedTFiles = [TFileNames[ind] for ind in sortedInds]
sortedTStrings = [TStrings[ind] for ind in sortedInds]
tStart = datetime.now()

# Function to terminate the subprocess
def terminate_process(proc):
    try:
        # Terminate the whole process group
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        print(f"Terminated subprocess with PID {proc.pid}")
    except Exception as e:
        print(f"Error terminating subprocess with PID {proc.pid}: {e}")


png_file_vec1=[one_dir+"/avg_polarization.png" for one_dir in sortedTFiles]
png_file_vec2=[one_dir+"/dipole_each_site.png" for one_dir in sortedTFiles]

for file_path in png_file_vec1:
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
            print(f"File {file_path} has been removed successfully.")
        except Exception as e:
            print(f"Error occurred while removing the file: {e}")
    else:
        print(f"File {file_path} does not exist.")

for file_path in png_file_vec2:
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
            print(f"File {file_path} has been removed successfully.")
        except Exception as e:
            print(f"Error occurred while removing the file: {e}")
    else:
        print(f"File {file_path} does not exist.")

