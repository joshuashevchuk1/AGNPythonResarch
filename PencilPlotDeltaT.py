import pencil_old as pc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys
import traceback

CONST_INTERVAL=100
max_orbits = 500
inital_ivar = 0

def plots():
    root = os.getcwd()  # root dir is fucked up due to a space
    dir_run_list = next(os.walk('.'))[1]
    for i in range(len(dir_run_list)):
        print("dir_run_list is " + dir_run_list[i])
        os.chdir(dir_run_list[i])
        print("current cwd is " + os.path.split(os.getcwd())[1])
        try:
            plotRuns(inital_ivar)
            os.chdir(root)
        except:
            print("============")
            print("ignoring run")
            print("============")
            traceback.print_exc()
            os.chdir(root)
    os.chdir(root)

def plotRuns(ivar):
    print("============")
    print("entering plotRuns")
    DTarray=[]
    while ivar <= max_orbits:
        try:
            DTarray = getRunData(ivar,DTarray)
        except:
            print("bad run or no run")
            traceback.print_exc()
        ivar = ivar + CONST_INTERVAL
    plotCollectedData(DTarray)
    print("leaving plotRuns")
    print("============")

def getRunData(ivar, paramDTarray):
    ff = pc.read_var(trimall=True, ivar=ivar, magic=["TT"],quiet=True)
    ff0 = pc.read_var(trimall=True, ivar=0, magic=["TT"],quiet=True)
    dfT = ff.TT[:] - ff0.TT[:]
    paramDTarray.append(np.mean(np.log(np.sum(dfT ** 2, axis=0))))
    return paramDTarray

def plotCollectedData(paramDTarray):
    print("entering plotCollectedData")
    plt.plot(paramDTarray)
    plt.grid(True)
    dir_gamma_patches = mpatches.Patch(
        color='white',
        label=r'$\gamma$ :'
              + str(paramDTarray))
    plt.legend(handles=[dir_gamma_patches], loc=2)
    plt.savefig("DeltaT" + "-temp-" + ".png")
    plt.close()
    print("leaving plotCollectedData")

plots()
