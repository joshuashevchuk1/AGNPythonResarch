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

while wait <= 2:
    while i <= len(dir_run_list)-1:
        try:
            os.chdir(str(dir_run_list[i]))
            # print('===================')
            #print('changing directory')
            # print('===================')
            os.system('sbatch STest')
            os.chdir('..')
        except:
            print('===================')
            print('directory not found')
            print('===================')
            traceback.print_exc()
        i = i+di
    i = 0
    print('===================')
    print('waiting for jobs to finish')
    print('===================')
    time.sleep(5)
    os.system('squeue -u joshshev')
    time.sleep(300)
    print('===================')
    print('adding time delay for nodes to stop')
    print('===================')
    time.sleep(60)
    os.system('squeue -u joshshev')
    time.sleep(5)
    # time.sleep(5)
    wait = wait+dwait
