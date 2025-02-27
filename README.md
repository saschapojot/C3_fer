In the very beginning, initialize directories:
python mk_dir.py, to set coefficients, T, and directories
##########################################
To manually perform each step of computations
1. python launch_one_run.py ./path/to/mc.conf
2. make run_mc
3. ./run_mc ./path/to/cppIn.txt
4. python check_after_one_run.py ./path/to/mc.conf lastFileNum
5. go to 1, until no more data points are needed

#########################################
To run 1 pass of mc with checking statistics
1. cmake .
2. make run_mc
3. python exec_checking.py T N lastFileNum