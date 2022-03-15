import os
import sys

dir_run_list = next(os.walk('.'))[1]

print(dir_run_list)

i = 0
di = 1

wait = 0
dwait = 1

while i <= len(dir_run_list)-1:
    os.system('cp -rf SContinue '+str(dir_run_list[i]))
    os.system('cp -rf SBatch '+str(dir_run_list[i]))
    os.system('cp -rf SContinue_2hrs '+str(dir_run_list[i]))
    os.system('cp -rf SBatch_2hrs '+str(dir_run_list[i]))
    os.system('cp -rf SDev '+str(dir_run_list[i]))
    os.system('cp -rf SBatch_5 '+str(dir_run_list[i]))
    os.system('cp -rf SContinue_24hrs '+str(dir_run_list[i]))
    os.system('cp -rf *.py '+str(dir_run_list[i]))
    i = i+di
