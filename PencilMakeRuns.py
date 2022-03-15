#
# test out dir making script
#

import sys
import os

dir_run_list = next(os.walk('.'))[1]

print(dir_run_list)
print(dir_run_list[0])

# 3 is a placeholder


length = 2
run_start = 900
run_step = 1
run_length = int(run_start)

NewDirs = []


os.chdir(dir_run_list[1])

while run_start <= length+run_length:
    os.system('pc_newrun AGNRun'+str(run_start-1+run_step))
    NewDirs.append(run_start)
    run_start = run_start+run_step

os.chdir('..')

print(NewDirs)
