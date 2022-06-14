import pencil_old as pc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys
import traceback

def plots():
    root = os.getcwd()  # root dir is fucked up due to a space
    dir_run_list = next(os.walk('.'))[1]
    for i in range(len(dir_run_list)):
        print("dir_run_list is " + dir_run_list[i])
        os.chdir(dir_run_list[i])
        print("current cwd is " + os.path.split(os.getcwd())[1])
        try:
            plotRuns(str(os.path.split(os.getcwd())[1]))
            os.chdir(root)
        except:
            print("============")
            print("ignoring run")
            print("============")
            traceback.print_exc()
            os.chdir(root)
    os.chdir(root)


def plotRuns(dirName):
    print("============")
    print("entering plotRuns")
    plotCollectedData(getScaleRatio(),dirName)
    print("leaving plotRuns")
    print("============")

def getScaleHeight():
    print('entering getAspectRatio')
    ts = pc.read_ts()
    t = ts.t

    radius = ts.xq2

    par = pc.read_param()
    gravC = par.g0
    Mstar = par.pmass[0]
    cs = par.cs0
    gamma = par.gamma

    Kepler_F = np.sqrt(gravC * Mstar / radius)
    aspect_ratio = cs / Kepler_F
    scale_height = aspect_ratio*radius
    scale_height = scale_height[:getCutOff()]
    print('Leaving getAspectRatio')
    return scale_height

def getCutOff():
    global q
    global ecc_int
    print('entering getCutOff')
    ts = pc.read_ts()
    t = ts.t

    radius = ts.xq2
    LinearVelocity = ts.vxq2
    AngularVelocity = ts.vyq2

    par = pc.read_param()

    if (par.iprimary == 1):
        q = par.pmass[1]
    else:
        q = par.pmass[0]

    v2 = LinearVelocity ** 2 + AngularVelocity ** 2
    semi_major = 1. / (2 / radius - v2)
    DArclength = radius ** 2 * (AngularVelocity / radius)
    ep1 = (DArclength ** 2) / semi_major
    eccentricity = (1 - ep1) ** 0.5
    ecc_int = par.eccentricity
    ecc = eccentricity

    indexTimeCutOff = 0

    for i in range(len(ecc)):
        if ecc_int != 0:
            if ecc[i] <= 0.01:
                indexTimeCutOff = i
                break

    timeCutOff = t[indexTimeCutOff] / (np.pi * 2)
    print('timeCutOff is ',np.round(timeCutOff))
    print('leaving get cutoff')
    return indexTimeCutOff

def plotCollectedData(paramDTarray, dir):
    print("entering plotCollectedData")
    print("paramDTArray is ",paramDTarray)
    plt.plot(paramDTarray)
    plt.grid(True)
    Dir_max_patches = mpatches.Patch(
        color='white', label=r'$h_{max}$ :' + str(max(paramDTarray)))
    Dir_min_patches = mpatches.Patch(
        color='white', label=r'$h_{min}$ :' + str(min(paramDTarray)))
    plt.legend(handles=[Dir_max_patches,Dir_min_patches], loc=2)
    plt.title('q = ' + str(q) + ',' + r'$\varepsilon$ = ' + str(ecc_int))
    plt.xlabel(r'$t/T_0$')
    plt.ylabel('h')
    plt.tight_layout()
    plt.savefig("Normal-h-new-" + str(dir) + "-temp-max-" + ".png")
    plt.close()
    print("leaving plotCollectedData")

def run():
    plots()

CONST_INTERVAL = 15
inital_ivar = 1
run()
