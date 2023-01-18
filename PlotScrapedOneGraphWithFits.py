import pencil_old as pc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys
import traceback
import json

scrapeDict={}
max_orbits = int(input('max Orbits? : '))


def plotCollectedData(paramDTarray):
    print("entering plotCollectedData")
    print("paramDTArray is ",paramDTarray)
    plt.plot(paramDTarray)
    print("leaving plotCollectedData")

def setScrapedJson():
    global scrapeDict
    with open('ScrapeDeltaTData_At'+str(max_orbits)+'.json') as f:
        scrapeDict = json.load(f)

def plot():
    global scrapeDict

    for key in scrapeDict:

        N=50

        kernel = np.ones((N,))/N

        temp = np.gradient(scrapeDict[key]["DTarray"])
        temp = np.convolve(temp,kernel,mode='valid')
            

        qnum = scrapeDict[key]["q"]

        q = "; q = " + str(scrapeDict[key]["q"])
        ecc = r"$\varepsilon$ = " + str(scrapeDict[key]["ecc_int"])

        if (qnum == 2e-4):
            plt.plot(temp, label=ecc + q,ls=":")
        else: 
            plt.plot(temp, label=ecc + q)

    A = [3.2460232783611094,
         3.2504927293831987,
         2.7923767202168905,
         2.7303702980173794,
         2.755121883800762]

    B = [-13.887087783571568,
          -13.950474787263751,
          -12.71500782691691,
          -12.322767961425448,
          -12.427188684084623]

    time = np.arange(0,175,1)

    for i in range(len(A)):
        f = A[i] * np.log(time) + B[i]
        plt.plot(f,time,label="i = " + str(i))

    plt.grid(True)
    plt.legend()
    plt.title('data')
    plt.xlabel(r'$t/T_0$')
    plt.ylabel(r'$\nabla'+'(Temperature)')
    plt.ylim(0,0.05)
    plt.tight_layout()
    plt.savefig("ScrapedDeltaData_" + str(max_orbits) + "_plot" + ".png")
    plt.close()

    print("entering plot")

setScrapedJson()
plot()
