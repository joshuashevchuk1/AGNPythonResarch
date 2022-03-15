#
# test out dir making script
#

import sys
import os

dir_run_list = next(os.walk('.'))[1]

print(dir_run_list)
print(dir_run_list[1])

Run_to_copy = dir_run_list[1]

i = 0
di = 1

while i <= len(dir_run_list)-1:
    os.chdir(str(dir_run_list[i]))
    os.system('pc_build -H login1.stampede2.tacc.utexas.edu.conf')
    os.chdir('..')
    i = i+di
