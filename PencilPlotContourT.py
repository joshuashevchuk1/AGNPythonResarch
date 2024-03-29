import pencil_old as pc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys
import traceback
import json

lastPointArray=[]
eccIntArray=[]
qArray=[]
q=0
ecc_int=0

def plots():
    root = os.getcwd()  # root dir is fucked up due to a space
    dir_run_list = next(os.walk('.'))[1]
    for i in range(len(dir_run_list)):
        print("dir_run_list is " + dir_run_list[i])
        os.chdir(dir_run_list[i])
        print("current cwd is " + os.path.split(os.getcwd())[1])
        try:
            addLastPoint(inital_ivar)
            os.chdir(root)
        except:
            print("============")
            print("ignoring run ",os.path.split(os.getcwd())[1])
            print("============")
            traceback.print_exc()
            os.chdir(root)
    os.chdir(root)

    print("lastPointArray is ",lastPointArray)
    print("eccIntArray is ",eccIntArray)
    print("qArray is ",qArray)

    plotContourRates()

def plotContourRates():
    global lastPointArray
    global qArray
    global eccIntArray
    name = os.path.split(os.getcwd())[1]
    ncolors=256

    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
    fig.subplots_adjust(bottom=0.07, top=0.95)
    PL2 = ax1.contour3D(eccIntArray, qArray, lastPointArray, ncolors)
    ax1.set_aspect('equal')
    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    cax.set_aspect(20)
    cax.set_ylabel('dT/dt', fontsize=10)
    plt.colorbar(PL2, cax=cax)
    plt.xlabel(r'$\epsilon')
    plt.ylabel('q')

    plt.savefig(name + "-RateContour-" + ".png")

def addLastPoint(ivar):
    global lastPointArray
    global qArray
    global eccIntArray
    print("============")
    print("entering plotRuns")

    timeCutOff = getCutOff()
    max_orbits = np.round(timeCutOff)

    DTarray = []
    while ivar <= max_orbits:
        try:
            DTarray = getRunData(ivar, DTarray)
        except:
            print("bad run or no run")
            traceback.print_exc()
        ivar = ivar + CONST_INTERVAL

    lastPoint = getRateOfLastPoint(DTarray)
    lastPointArray.append(lastPoint)
    qArray.append(q)
    eccIntArray.append(ecc_int)

    print("leaving plotRuns")
    print("============")


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
    return timeCutOff

def getRunData(ivar, paramDTarray):
    ff = pc.read_var(trimall=True, ivar=ivar, magic=["TT"], quiet=True)
    ff0 = pc.read_var(trimall=True, ivar=0, magic=["TT"], quiet=True)
    dfT = ff.TT[:] - ff0.TT[:]

    paramDTarray.append(np.max(np.log(np.sum(dfT ** 2, axis=0))))
    return paramDTarray

def getRateOfLastPoint(paramDTarray):
    gradientArray = np.gradient(paramDTarray)
    lastRate = gradientArray[len(gradientArray)-1]
    return lastRate


CONST_INTERVAL = 5
inital_ivar = 1
plots()
