import numpy
import sys
import traceback
import os
import logging
import matplotlib.pyplot

# run analysis for all AGNRuns

os.chdir('./Analysis')
data = os.listdir('.')

print(data)

for AGNRun in range(len(data)):
    print(str(data[AGNRun]))
    DIR = str(data[AGNRun])
    os.chdir(DIR)

    print('=========')
    print('=========')
    print('Analyzing next run')
    print('=========')
    print('=========')

    os.system('ipython EccDecayCN2008_editWL.py')
    os.system('ipython Torque.py')
    os.system('git add *.png')
    os.system('git add *.in')
    os.system('git add src/cparam.local')
    os.system('git add src/Makefile.local')
    os.system('git add *.py')
    os.chdir('..')

os.system('git commit')
os.system('git push')
