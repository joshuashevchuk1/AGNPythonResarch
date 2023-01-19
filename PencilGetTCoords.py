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
scrapeDict={}
max_orbits = int(input('max Orbits? : '))
timeCutOff=0

def scrape():
    root = os.getcwd()  # root dir is fucked up due to a space
    dir_run_list = next(os.walk('.'))[1]
    for i in range(len(dir_run_list)):
        print("dir_run_list is " + dir_run_list[i])
        os.chdir(dir_run_list[i])
        print("current cwd is " + os.path.split(os.getcwd())[1])
        try:
            addLastPoint(inital_ivar,str(os.path.split(os.getcwd())[1]))
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

    saveData()

def saveData():
    global scrapeDict
    data = scrapeDict
    with open('PencilGetTCoords'+str(max_orbits)+'.json', 'w') as f:
        json.dump(data, f)

def getCutOff():
    global q
    global ecc_int
    global timeCutOff
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

def addLastPoint(ivar,dir):
    global lastPointArray
    global qArray
    global eccIntArray
    global scrapeDict
    global max_orbits
    global timeCutOff

    print("============")
    print("entering plotRuns")

    initPars()

    DTarray = []
    try:
        DTarray = getRunData(ivar, DTarray)
    except:
        print("bad run or no run")
        traceback.print_exc()

    qArray.append(q)
    eccIntArray.append(ecc_int)
    localDict={"q":q,"ecc_int":ecc_int,"timeCutOff":timeCutOff,"DTarray":DTarray}

    scrapeDict[str(dir)] = localDict

    print("leaving plotRuns")
    print("============")
def initPars():
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
    ff1 = pc.read_var(trimall=True, ivar=1, magic=["TT"], quiet=True)
    ff0 = pc.read_var(trimall=True, ivar=0, magic=["TT"], quiet=True)
    dfT = ff.TT[:] - ff0.TT[:]
    df1 = ff.TT[:] - ff0.TT[:]

    paramDTarray.append(np.log(np.sum(df1 ** 2, axis=0)))
    paramDTarray.append(np.log(np.sum(dfT ** 2, axis=0)))
    plt.plot(paramDTarray[0])
    plt.plot(paramDTarray[1])
    plt.show()
    return paramDTarray


CONST_INTERVAL = 1
inital_ivar = 1
scrape()