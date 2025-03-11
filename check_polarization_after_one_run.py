import re
import subprocess
import sys
from datetime import datetime
import json
argErrCode=2

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    print("example: python check_polarization_after_one_run.py ./path/to/mc.conf lastFileNum")
    exit(argErrCode)

confFileName=str(sys.argv[1])
# print("confFileName is "+confFileName)
lastFileNum=int(sys.argv[2])

print(f"lastFileNum={lastFileNum}")
invalidValueErrCode=1
summaryErrCode=2
loadErrCode=3
confErrCode=4

#################################################
#parse conf, get jsonDataFromConf
confResult=subprocess.run(["python3", "./init_run_scripts/parseConf.py", confFileName], capture_output=True, text=True)
confJsonStr2stdout=confResult.stdout
# print(confJsonStr2stdout)
if confResult.returncode !=0:
    print("Error running parseConf.py with code "+str(confResult.returncode))
    # print(confResult.stderr)
    exit(confErrCode)


match_confJson=re.match(r"jsonDataFromConf=(.+)$",confJsonStr2stdout)
if match_confJson:
    jsonDataFromConf=json.loads(match_confJson.group(1))
else:
    print("jsonDataFromConf missing.")
    exit(confErrCode)
# print(jsonDataFromConf)
################################################

##########################################################
#statistics
t_stats_start=datetime.now()
checkU_dipole_ErrCode = 5
checkU_dipole_Process = subprocess.Popen(
    ["python3", "-u", "./oneTCheckObservables/check_dipole_OneT_pkl.py"
        , json.dumps(jsonDataFromConf),str(lastFileNum)],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
)
# Read output in real-time
while True:
    output = checkU_dipole_Process.stdout.readline()
    if output == '' and checkU_dipole_Process.poll() is not None:
        break
    if output:
        print(output.strip())

# Collect remaining output and error messages
stdout, stderr = checkU_dipole_Process.communicate()
# Check if the process was killed
if checkU_dipole_Process.returncode is not None:
    if checkU_dipole_Process.returncode < 0:
        # Process was killed by a signal
        print(f"checkU_dipole_Process was killed by signal: {-checkU_dipole_Process.returncode}")
    else:
        # Process exited normally
        print(f"checkU_dipole_Process exited with return code: {checkU_dipole_Process.returncode}")
else:
    print("checkU_dipole_Process is still running")
# Print any remaining standard output
if stdout:
    print(stdout.strip())

# Handle errors and print the return code if there was an error
if stderr:
    print(f"checkU_dipole_Process return code={checkU_dipole_Process.returncode}")
    print(stderr.strip())
t_stats_End=datetime.now()

print(f"stats time: {t_stats_End-t_stats_start}")
##########################################################