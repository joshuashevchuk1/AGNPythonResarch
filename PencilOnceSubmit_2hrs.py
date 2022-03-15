import numpy as np
import os
import sys
import traceback
import time

# submit jobs on stampede2 by using every directorys batch file

print(os.getcwd())

dir_run_list = next(os.walk('.'))[1]

print(dir_run_list)

i = 0
di = 1

wait = 0
dwait = 1

while i <= len(dir_run_list)-1:
    try:
        os.chdir(str(dir_run_list[i]))
        # print('===================')
        #print('changing directory')
        # print('===================')
        try:
            os.system('sbatch SBatch_5')
        except:
            traceback.print_exc()
            os.chdir('..')
        os.chdir('..')
    except:
        print('===================')
        print('directory not found')
        print('===================')
        traceback.print_exc()
    i = i+di
print('===================')
print('Submitting Batch Jobs')
print('===================')
os.system('squeue -u joshshev')
print('===================')
print('adding time delay for nodes to stop')
print('===================')
os.system('squeue -u joshshev')
