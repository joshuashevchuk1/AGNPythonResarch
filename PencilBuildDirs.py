#
# test out dir making script
#

import sys
import os

dir_run_list = next(os.walk('.'))[1]

print(dir_run_list)
print(dir_run_list[1])

Run_to_copy = dir_run_list[1]

length = 2
run_start = 825
run_step = 1
run_length = int(run_start)

NewDirs = []

while run_start <= length+run_length:
    os.system('mkdir AGNRun'+str(run_start-1+run_step))
    print('Making Dir')
    NewDirs.append('AGNRun'+str(run_start))
    run_start = run_start+run_step

i = 0
di = 1

while i <= len(NewDirs)-1:
    os.chdir(str(NewDirs[i]))
    os.system('pc_setupsrc')
    os.chdir('..')
    i = i+di
i = 0
di = 1

while i <= len(NewDirs)-1:
    os.system('cp -rf '+str(Run_to_copy)+'/*.in '+str(NewDirs[i]))
    os.system('cp -rf '+str(Run_to_copy) +
              '/src/cparam.local '+str(NewDirs[i])+'/src/')
    os.system('cp -rf '+str(Run_to_copy) +
              '/src/Makefile.local '+str(NewDirs[i])+'/src/')
    i = i+di

i = 0
di = 1

while i <= len(NewDirs)-1:
    os.chdir(str(NewDirs[i]))
    os.system('mkdir data')
    os.chdir('..')
    i = i+di

i = 0
di = 1


while i <= len(NewDirs)-1:
    os.chdir(str(NewDirs[i]))
    os.system('pc_build -H login1.stampede2.tacc.utexas.edu.conf')
    os.chdir('..')
    i = i+di
